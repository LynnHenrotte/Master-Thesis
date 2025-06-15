import sys
import statistics

quality_scores = sys.argv[1]

with open(quality_scores, "r") as qual:
    # Extract quality scores without newlines
    scores = qual.read().splitlines()
    
    if len(scores) > 0:
        # Convert strings to floats
        scores = [float(score) for score in scores]

        # Compute mean quality score and print to standard output
        mean_score = round(statistics.mean(scores), 2)
        print(mean_score)
    
    else:
        print(0)
