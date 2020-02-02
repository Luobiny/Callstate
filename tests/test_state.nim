import ../callable

var raw_depth: seq[int32] =      @[0i32, 30i32, 30i32, 16i32, 8i32, 30i32, 30i32, 40i32, 5000i32, 5i32]
var qc_depth: seq[int32] =       @[0i32, 22i32, 19i32, 6i32,  2i32, 20i32, 21i32, 2i32, 4000i32, 2i32]
var low_mapq_depth: seq[int32] = @[0i32, 8i32,  10i32, 10i32, 6i32, 9i32, 8i32, 30i32, 10i32, 1i32]

proc test_get_state() =
  for i, raw in raw_depth:
    var qc = qc_depth[i]
    var low_mapq = low_mapq_depth[i]
    var state = get_state(raw, qc, low_mapq, min_depth=20, min_mapq=10, max_depth=2000,
                          min_depth_low_mapq=10, max_low_mapq_frac=0.5)
    case i:
        of 0:
            assert state == CalledStates.NO_COVERAGE
        of 1, 5, 6:
            assert state == CalledStates.CALLABLE
        of 2, 4, 9:
            assert state == CalledStates.LOW_COVERAGE
        of 3, 7:
            assert state == CalledStates.POOR_MAPPING_QUALITY
        of 8:
            assert state == CalledStates.HIGH_COVERAGE
        else:
            assert true

proc test_gen_state() =
    let ref_N_bed = "ref_N.bed"
    var cv = (raw: raw_depth, qc: qc_depth, low_mapq: low_mapq_depth) 

    var regions = new_seq[region_t]()
    regions.add(region_t(chrom: "22", start: 0'u32, stop: 10'u32))

    var idx = 0
    for p in gen_state(cv, "22", regions, ref_N_bed):
        #echo p.start, ", ", p.stop, ", ", p.value
        case idx:
            of 0:
                assert (p.start, p.stop, p.value) == (0, 1, "NO_COVERAGE")
            of 1:
                assert (p.start, p.stop, p.value) == (1, 2, "CALLABLE")
            of 2:
                assert (p.start, p.stop, p.value) == (2, 3, "LOW_COVERAGE")
            of 3:
                assert (p.start, p.stop, p.value) == (3, 4, "POOR_MAPPING_QUALITY")
            of 4:
                assert (p.start, p.stop, p.value) == (4, 5, "LOW_COVERAGE")
            of 5:
                assert (p.start, p.stop, p.value) == (5, 7, "CALLABLE")
            of 6:
                assert (p.start, p.stop, p.value) == (7, 8, "POOR_MAPPING_QUALITY")
            of 7:
                assert (p.start, p.stop, p.value) == (8, 9, "HIGH_COVERAGE")
            of 8:
                assert (p.start, p.stop, p.value) == (9, 10, "REF_N")
            else:
                assert false
        idx += 1
 
    regions = new_seq[region_t]()
    regions.add(region_t(chrom: "22", start: 1'u32, stop: 3'u32))
    regions.add(region_t(chrom: "22", start: 6'u32, stop: 8'u32))
    idx = 0
    for p in gen_state(cv, "22", regions, ref_N_bed):
      echo p.start, ", ", p.stop, ", ", p.value
      case idx:
        of 0:
          assert (p.start, p.stop, p.value) == (1, 2, "CALLABLE")
        of 1:
          assert (p.start, p.stop, p.value) == (2, 3, "LOW_COVERAGE")
        of 2:
          assert (p.start, p.stop, p.value) == (6, 7, "CALLABLE")
        of 3:
          assert (p.start, p.stop, p.value) == (7, 8, "POOR_MAPPING_QUALITY")
        else:
          assert false 
      idx += 1


test_get_state()
test_gen_state()

