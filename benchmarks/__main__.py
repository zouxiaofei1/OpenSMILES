"""Allow `python -m benchmarks` as well as `python -m benchmarks.benchmark`."""

from benchmarks.benchmark import main

if __name__ == "__main__":
    main()
