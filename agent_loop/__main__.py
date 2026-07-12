"""Allow `python -m agent_loop` as well as `python -m agent_loop.loop`."""

from agent_loop.loop import main

if __name__ == "__main__":
    main()
