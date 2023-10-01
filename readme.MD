# Genome to transcript mapper
Convert coordinates or intervals from genome space to transcript space.

## Installation
Via conda: `conda install -c molinerislab genometotranscriptmapper`

Via pip: clone the repository, run `pip install . `

# Usage inside python scripts
GenomeToTranscriptMapper class can be imported inside python scripts: from genomeToTranscriptMapper import GenomeToTranscriptMapper
## __init__:
Initialize mapper object with the transcript informations:
* exons_interval: list of tuples in the form (begin, end) representing exons
* * intervals are left inclusive and right exclusive.
* * order is not important
* strand: "+" or "-"
    
## convert_points_genome_to_transcript:
Convert a list of points from genome space to transcript space
Please note: if you have many coordinates to convert, calling the method once with a list of points is way faster than calling it many times with single values.
Parameters:
* positions: list of integers representing genome space positions.
* outside_exons_points [optional]: allow to specify the behaviour for points outside of transcript:
* * OutOfExonsPoint.THROW_EXCEPTION [default behaviour]
* * OutOfExonsPoint.RETURN_MINUS_ONE
* * OutOfExonsPoint.RESTRICT_SMALLER: return first transcript position on the left of the specified one.
* * OutOfExonsPoint.RESTRICT_BIGGER: return first transcript position on the right of the specified one.

## convert_interval_genome_to_transcript:
Convert an interval from genome space to transcript space.  
Parameters:
* begin, end: integers.
* restrict_interval_to_transcript: specify what to do if the interval is only partially
*   included in the transcript:
* * True: interval is cut to fit inside the transcript.
* * False: if interval isn't completely included in the transcript, an execption is raised.
NOTE: all intervals are considered left side inclusive and right side exclusive
For example, for strand + : (1000, 2000) => from 1000 to 1999 included
for strand - : (2000, 1000) => from 2000 to 1001 included
or, strand - reversed: (1000, 2000) => from 1999 to 1000


# Usage from command line
The tool can be used from command line:
* stdin: intervals to convert, in bed format
* argument: exon intervals in bed format
Ex: `cat test/to_convert_example.bed | genomeToTranscriptMapper test/transcript_exons_example.bed`

Note regarding the bed format:
genomeToTranscriptMapper need the first 6 fields of bed files; Other fields will be preserved in the output.
