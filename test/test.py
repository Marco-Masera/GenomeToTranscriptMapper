import sys
import os
current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)
from src import genomeToTranscriptMapper

def run_tests():
    #Positive strand tests
    transcript = genomeToTranscriptMapper.GenomeToTranscriptMapper([(0, 101), (200, 301)], "+")
    assert transcript.convert_points_genome_to_transcript([10]) == 10
    assert transcript.convert_points_genome_to_transcript([210]) == 111
    assert transcript.convert_points_genome_to_transcript([150], genomeToTranscriptMapper.OutOfExonsPoint.RETURN_MINUS_ONE) == -1 
    assert transcript.convert_points_genome_to_transcript([150], genomeToTranscriptMapper.OutOfExonsPoint.RESTRICT_SMALLER) == 100
    assert transcript.convert_points_genome_to_transcript([150], genomeToTranscriptMapper.OutOfExonsPoint.RESTRICT_BIGGER) == 101
    assert transcript.convert_interval_genome_to_transcript(200, 202) == [101, 103]

    #Negative strand tests
    transcript2 = genomeToTranscriptMapper.GenomeToTranscriptMapper([(30, 19), (10, -1)], "-")
    assert transcript2.convert_points_genome_to_transcript([30]) == 0
    assert transcript2.convert_points_genome_to_transcript([20]) == 10
    assert transcript2.convert_points_genome_to_transcript([10]) == 11
    assert transcript2.convert_interval_genome_to_transcript(10, 21) == [10, 12]
    
if __name__=="__main__":
    run_tests()