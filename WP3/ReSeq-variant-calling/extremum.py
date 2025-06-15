import sys

select = sys.argv[1]
numbers_file = sys.argv[2]  # File containing a single (floating point) number on each line

try:
    with open(numbers_file, "r") as file:

        # Extract numbers from file
        numbers = [float(line.rstrip()) for line in file]

        # Compute maximum or minimum of numbers
        if select == "max":
            print(max(numbers))
        elif select == "min":
            print(min(numbers))
        else:
            print("Select one of the following:\nmax: compute maximum\nmin: compute minimum")
except FileNotFoundError:
    print(0)
