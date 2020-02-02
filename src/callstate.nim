import hts
import docopt
import tables
import logging
import strutils
import strformat
import ./refn

var logger = newConsoleLogger(fmtStr="[$time] - $levelname: ",
                              levelThreshold=lvlInfo)

type
  item = tuple[start: int, stop: int, offset: int]
  depth_s = tuple[start: int, stop: int, value: string]

  coverage_t* = seq[int32]
  cov_t* = tuple[raw: coverage_t, qc: coverage_t, low_mapq: coverage_t]

  CalledStates* = enum
    NONE = -1, REF_N = 0, NO_COVERAGE = 1, LOW_COVERAGE = 2, CALLABLE = 3, 
    EXCESSIVE_COVERAGE = 4, POOR_MAPPING_QUALITY = 5


proc get_state*(raw: int, qc: int, low_mapq: int, min_depth: int,
               min_mapq: int, max_depth: int, min_depth_low_mapq: int,
               low_mapq_frac: float): CalledStates =
    # Given the raw depth, qc depth, and low_mapq depth of a locus, determine
    # the called states.
    if raw == 0:
      return CalledStates.NO_COVERAGE
    elif raw >= min_depth_low_mapq and (low_mapq / raw) >= low_mapq_frac :
      return CalledStates.POOR_MAPPING_QUALITY
    elif qc < min_depth:
      return CalledStates.LOW_COVERAGE
    elif raw >= max_depth and max_depth != -1:
      return CalledStates.EXCESSIVE_COVERAGE
    else:
      return CalledStates.CALLABLE


proc init(arr: var coverage_t, tlen: int) =
  # Try to re-use the array.
  if len(arr) != tlen:
    if arr.len == 0:
      arr = new_seq[int32](tlen)
      return
    else:
      arr.set_len(int(tlen))

  zeroMem(arr[0].addr, len(arr) * sizeof(arr[0]))
  shallow(arr)


iterator gen_state*(arr: cov_t, chrom: string, regions: seq[region_t],
                    min_depth: int=20, min_mapq: int=10, max_depth: int= -1,
                    min_depth_low_mapq: int=10,
                    low_mapq_frac: float=0.5): depth_s =
  # A generator that generates regions of different called state in the 
  # given regions.
  for region in regions:
    var
      state: CalledStates
      last_state = CalledStates.NONE
      last_ipos = -1

    for ipos in int(region.start) .. int(region.stop) - 1:
      if base_is_N(chrom, ipos):
          state = CalledStates.REF_N
      else:
        let raw: int = arr.raw[ipos]
        let qc: int = arr.qc[ipos]
        let low_mapq: int = arr.low_mapq[ipos]
        state = get_state(raw=raw, qc=qc, low_mapq=low_mapq, min_depth=min_depth,
                          min_mapq=min_mapq, max_depth=max_depth,
                          min_depth_low_mapq=min_depth_low_mapq,
                          low_mapq_frac=low_mapq_frac)
      
      # The last position
      if ipos == int(region.stop) - 1:
          if state == last_state:
              yield (last_ipos, ipos + 1, $last_state)
          else:
              yield (last_ipos, ipos, $last_state)
              yield (ipos, ipos + 1, $state) 
      else:
          if state == last_state:
              continue
          if last_state != CalledStates.NONE:
              yield (last_ipos, ipos, $last_state) 
     
          last_state = state
          last_ipos = ipos


proc get_tid(tgts: seq[hts.Target], chrom: string): int =
  for t in tgts:
    if t.name == chrom:
      return t.tid


iterator regions(bam: hts.Bam, region: region_t, tid: int, 
                 targets: seq[hts.Target]): Record {.inline.} =
  if region == nil:
    for r in bam:
      yield r
  elif region != nil:
    var stop = region.stop
    if tid != -1:
      if stop == 0:
        stop = targets[tid].length
      for r in bam.query(tid, int(region.start), int(stop)):
        yield r
    else:
      stderr.write_line(fmt"Error: {region.chrom} not found")


iterator gen_start_ends*(c: Cigar, ipos: int): item {.inline.} =
  # generate start, end, and offset.
  if c.len == 1 and c[0].op == CigarOp.match:
    yield (ipos, ipos + c[0].len, 0)
  else:
    var offset = 0
    var pos = ipos
    var con: Consume
    for op in c:
      con = op.consumes
      if not con.reference:
        if op.op == CigarOp.hardclip:
            continue
        else:
            offset += op.len
            continue
      var olen = op.len
      if con.query:
          yield (pos, pos + olen, offset)
      else:
        offset -= op.len
      pos += olen


proc inc_coverage(rec: Record, raw: var seq[int32], qc: var seq[int32], 
                  low_mapq: var seq[int32], min_mapq: int, max_low_mapq: int, 
                  min_base_qual: uint8) {.inline.} =

  var qual = new_seq[uint8]()

  for p in gen_start_ends(rec.cigar, rec.start.int):
      raw[p.start] += 1
      raw[p.stop]  -= 1
      if int(rec.mapping_quality) <= max_low_mapq:
          low_mapq[p.start] += 1
          low_mapq[p.stop]  -= 1
      if int(rec.mapping_quality) >= min_mapq:
          qc[p.start] += 1
          qc[p.stop]  -= 1
          if min_base_qual > 0'u8:
            qual = rec.base_qualities(qual)
            var s = p.start - rec.start + p.offset
            var e = p.stop  - rec.start + p.offset
            for idx, q in qual[s..<e]:
              if q < min_base_qual:
                qc[p.start + idx]     -= 1
                qc[p.start + idx + 1] += 1 


proc coverage(bam: Bam, arr: var cov_t, region: region_t, eflag: uint16=1796, 
              min_mapq: int=10, max_low_mapq: int=1,
              min_base_qual: uint8=10): int =
  var
    targets = bam.hdr.targets
    tgt: hts.Target

  var tid = if region != nil: get_tid(targets, region.chrom) else: -1
  if tid == -1:
    return -1

  tgt = targets[tid]
  var 
    found = false

  for rec in bam.regions(region, tid, targets):
    if not found:
      arr.raw.init(int(tgt.length+1))
      arr.qc.init(int(tgt.length+1))
      arr.low_mapq.init(int(tgt.length+1))
      found = true
    
    if (rec.flag and eflag) != 0:
      continue

    inc_coverage(rec, arr.raw, arr.qc, arr.low_mapq, min_mapq=min_mapq,
                 max_low_mapq=max_low_mapq, min_base_qual=min_base_qual)

  if not found:
    return -2
  return tgt.tid


proc to_coverage(c: var coverage_t) =
  # to_coverage converts from an array of start/end inc/decs to actual coverage.
    var cum = int32(0)
    for i, d in c.mpairs:
      cum += d
      d = cum

proc mean_depth*(depth: coverage_t): float =
    # Calculate the mean depth given the start and stop index of a region.
    var sum = 0
    for idx in 0 ..< len(depth):
        sum += depth[idx]
    return float(sum) / float(depth.len)


proc main(bam: hts.Bam, eflag: uint16, args: Table[string, docopt.Value]) =
  var
    arr: cov_t
    bed_regions: TableRef[string, seq[region_t]]
    r_chrom: region_t
    bedfile = $args["<BED>"]
    outfile = $args["--output"]
    targets = bam.hdr.targets
    min_base_qual      = parse_int($args["--min-base-qual"]).uint8
    min_mapq           = parse_int($args["--min-mapq"])
    min_depth          = parse_int($args["--min-depth"])
    max_depth          = parse_int($args["--max-depth"])
    max_low_mapq       = parse_int($args["--max-low-mapq"])
    low_mapq_frac      = parse_float($args["--low-mapq-frac"])
    min_depth_low_mapq = parse_int($args["--min-depth-low-mapq"])

  logger.log(lvlInfo, fmt"min-base-qual = {min_base_qual}")
  logger.log(lvlInfo, fmt"min-mapq = {min_mapq}")
  logger.log(lvlInfo, fmt"max-low-mapq = {max_low_mapq}")
  logger.log(lvlInfo, fmt"min-depth = {min_depth}")
  logger.log(lvlInfo, fmt"max-depth = {max_depth}")
  logger.log(lvlInfo, fmt"low-mapq-frac = {low_mapq_frac}")
  logger.log(lvlInfo, fmt"min-depth-low-mapq = {min_depth_low_mapq}")

  logger.log(lvlInfo, fmt"The BED file is: {bedfile}")

  bed_regions = bed_to_table(bedfile)
 
  var output: File
  if outfile != "" and outfile != "nil":
    logger.log(lvlInfo, fmt"The output file is: {outfile}")
    output = open($args["--output"], fmWrite)
  else: output = stdout
  let output_headers = @["#chrom", "start", "end",  "call_state", "raw_depth",
                         "qc_depth", "low_mapq_depth"]
  output.write_line(output_headers.join("\t"))

  let  bd_output = $args["--base-depth-output"]
  var bdout: File
  if bd_output != "" and bd_output != "nil":
    logger.log(lvlInfo, fmt"The per-base depth output file is: {bd_output}")  
    bdout = open(bd_output, fmWrite)
    bdout.write_line(output_headers.join("\t"))

  for target in targets:
    if not bed_regions.hasKey(target.name):
        logger.log(lvlDebug, fmt"Skip chrom {target.name}, it is not in {bedfile}")
        continue

    logger.log(lvlInfo, fmt"Generating coverage for chromosome {target.name}")
    r_chrom = region_t(chrom: target.name) 
    var tid = coverage(bam, arr, r_chrom, eflag=eflag,
                       min_base_qual=min_base_qual)
    if tid == -1: break # -1 means that chrom is not even in the bam
    if tid == -2: break # -2 means there were no reads in the bam
    logger.log(lvlInfo, "Calculating raw depth")
    arr.raw.to_coverage() 
    logger.log(lvlInfo, "Calculating QC depth")
    arr.qc.to_coverage() 
    logger.log(lvlInfo, "Calculating low MAPQ depth")
    arr.low_mapq.to_coverage() 

    logger.log(lvlInfo, "Writing callable states BED file...")
    for p in gen_state(arr, target.name, bed_regions[target.name], 
                       min_depth=min_depth, low_mapq_frac=low_mapq_frac):
      if bdout != nil:
        for ipos in p.start ..< p.stop:
          var bfs = @[target.name, $ipos, $(ipos + 1), p.value, $arr.raw[ipos],
                     $arr.qc[ipos], $arr.low_mapq[ipos]]
          bdout.write_line(bfs.join("\t"))

      # Write regular output file.      
      var fs = @[target.name, $p.start, $p.stop, p.value]
      let raw_d = mean_depth(arr.raw[p.start ..< p.stop])
      let qc_d = mean_depth(arr.qc[p.start ..< p.stop])
      let lmq_d = mean_depth(arr.low_mapq[p.start ..< p.stop])
      fs.add(fmt"{raw_d:.2f}")
      fs.add(fmt"{qc_d:.2f}")
      fs.add(fmt"{lmq_d:.2f}")
      output.write_line(fs.join("\t"))

  if output != stdout: output.close()
  if bdout != nil: bdout.close()
  logger.log(lvlInfo, "Done!")


when isMainModule:
  let version = "callstate 0.0.1"
  let doc = format("""
  $version

  Usage: callstate [options] <BED> <BAM>

Arguments:

  <BAM>  the alignment file for which to calculate callable states
  <BED>  The BED file that contains regions.

Common Options:

  -t --threads <threads>                 Number of BAM decompression threads [default: 4]
  -o --output <output>                   The output BED file
  -bdout --base-depth-output <bdout>     If a file name is given, per-base depth will be written to this file

Quality Metrics Options:

  -mbq --min-base-qual <mbq>             The minimum base quality for a base to contribute to QC depth [default: 10]
  -mmq --min-mapq <mmq>                  The minimum mapping quality of reads to count as QC depth [default: 10]
  -mdp --min-depth <mdp>                 The minimum QC read depth before a read is considered callable [default: 20]
  -mlmq --max-low-mapq <mlmq>            The maximum value of MAPQ before a read is considered as problematic mapped read. [default: 1]
  -mxdp --max-depth <mxdp>               The maximum read depth before a locus is considered high coverage [default: -1]
  -mdflmq --min-depth-low-mapq <mdflmq>  Minimum read depth before a locus is considered candidate for poorly mapped [default: 10]
  -frlmq --low-mapq-frac <frlmq>         If the fraction of low mapping reads exceeds this value, the site is considered poorly mapped [default: 0.5]

Other options:
  
  -F --flag <FLAG >                      exclude reads with any of the bits in FLAG set [default: 1796]
  -h --help                              show help
  """ % ["version", version])

  var args: Table[string, Value]
  try:
    args = docopt(doc, version = version, quit=false)
  except DocoptExit:
    echo (ref DocoptExit)(get_current_exception()).usage
    quit -1

  var 
    bam: Bam
    eflag: uint16 = uint16(parse_int($args["--flag"]))
    bamfile = $args["<BAM>"]
    threads = parse_int($args["--threads"])

  logger.log(lvlInfo, fmt"The alignment file is: {bamfile}")
  open(bam, bamfile, threads=threads, index=true)
  if bam.idx == nil:
    logger.log(lvlError, fmt"Alignment file {bamfile} must be indexed")
    quit(2)

  var opts = SamField.SAM_FLAG.int or SamField.SAM_RNAME.int or 
             SamField.SAM_POS.int or SamField.SAM_MAPQ.int or 
             SamField.SAM_CIGAR.int
  discard bam.set_option(FormatOption.CRAM_OPT_REQUIRED_FIELDS, opts)
  discard bam.set_option(FormatOption.CRAM_OPT_DECODE_MD, 0)

  main(bam=bam, eflag=eflag, args)

