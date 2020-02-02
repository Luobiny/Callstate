# Callstate
A tool that generates similar results as CallableLoci in GATK3 but runs much faster.

## Algorithm
This tool uses similar depth calculating algorithm in Mosdepth to calculate the raw depth, qc depth and low_mapq depth and then uses the same algorithm in CallableLoci tool in GATK3 to determine the following callable states:

 - REF_N
 - NO_COVERAGE
 - LOW_COVERAGE
 - POOR_MAPPING_QUALITY
 - EXCESSIVE_COVERAGE
 - CALLABLE

For each region in the output file, besides the location of region, Callstate also outputs the mean raw depth, mean QC depth, and the mean Low MAPQ depth 

## Usage
```
callstate 0.0.1

  Usage: callstate [options] <BED> <BAM>

Arguments:

  <BAM>  the alignment file for which to calculate callable states
  <BED>  The BED file that contains regions.

Common Options:

  -t --threads <threads>                 Number of BAM decompression threads [default: 4]
  -o --output <output>                   The output BED file.
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
```

## Examples
To run with default values for all parameters:
```
callstate -o <output.bed> <capture.bed> <bam/final.bam>
```

## Sample output
```
#chrom  start   end     call_state      raw_depth       qc_depth        low_mapq_depth
22      17565780        17566404        CALLABLE        309.22  305.61  0.00
22      17577753        17577787        LOW_COVERAGE    10.41   9.76    0.00
22      17577787        17577789        POOR_MAPPING_QUALITY    161.78  31.23   129.76
22      17577789        17577790        LOW_COVERAGE    20.00   19.00   0.00
22      17577790        17578173        CALLABLE        120.22  119.55  0.00
```
