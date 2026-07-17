import pytest
import logging
import openawsem.scripts.mm_run as mm_run
import math

# Define parameterized test cases
@pytest.mark.parametrize(
    ["steps", "numFrames", "reportInterval", "should_raise", "expected_steps", "expected_frames", "expected_interval"], 
    [
        # Test 1: Default, numFrames = 400, steps = 1e7 -> steps / numFrames
        (None, None, None, False, 1e7, 400, math.ceil(1e7 / 400)),
        # Test 2: steps = 100, interval = 99, adjust steps to 198 (ceil)
        ("100", None, "99", False, 198, 2, 99),
        # Test 3: steps = 100, interval = 100, 1 frame recorded
        ("100", None, "100", False, 100, 1, 100),
        # Test 4: interval > steps (101 > 100), should adjust steps to 101
        ("100", None, "101", False, 101, 1, 101),
        # Test 5: steps = 100, numFrames = 99, adjust steps to match numFrames
        ("100", "99", None, False, 198, 99, 2),
        # Test 6: steps = 100, numFrames = 100, interval should be 1
        ("100", "100", None, False, 100, 100, 1),
        # Test 7: numFrames > steps (101 > 100), adjust number of frames to match steps
        ("100", "101", None, False, 100, 100, 1),
        # Test 8: Conflicting numFrames and interval, Interval is prefered
        ("100", "100", "100", False, 100, 1, 100),
        # Test 11: Interval = 1, should work, 1 frame per step
        ("100", "100", "1", False, 100, 100, 1),
        # Test 12: Single frame at the end, interval = 100, 1 frame
        ("100", "1", "100", False, 100, 1, 100),
        # Test 13: steps = 500, numFrames = 400
        ("500", None, None, False, math.ceil(500 / 400)*400, 400, math.ceil(500 / 400)),
        # Test 14: steps = 500, numFrames = 101, adjust steps to 505 (ceil)
        ("500", "101", None, False, 505, 101, 5),
        # Test 15: steps = 500, interval = 101, adjust steps to 505 (ceil)
        ("500", None, "101", False, 505, 5, 101),
        # Test 16: Invalid number of frames, should raise error
        ("10", "-1", None, True, None, None, None),
        # Test 17: steps = 0, should adjust steps to 1e7
        ("0", None, None, False, 1e7, 400, math.ceil(1e7 / 400)),
        # Test 18: Invalid interval, should raise error
        ("-100", None, "-1", True, None, None, None),
    ]
)
def test_numFrames_reportFrequency_handling(
        steps, numFrames, reportInterval, should_raise, expected_steps, expected_frames, expected_interval):
    
    # Convert inputs to command line arguments for argparse
    cmd_args = [
        "test_script.py",  # Simulating script name
        "1r69",  # PDB ID (protein name)
        "--dryRun",  # Dry run option to return args without running the simulation
    ]
    
    # Add optional arguments based on test inputs
    if steps:
        cmd_args.extend(["--steps", steps])
    if numFrames:
        cmd_args.extend(["--numFrames", numFrames])
    if reportInterval:
        cmd_args.extend(["--reportInterval", reportInterval])

    print(f"Running with arguments: {cmd_args}")

    # Test whether ValueError should be raised or not
    if should_raise:
        with pytest.raises(ValueError):
            mm_run.main(cmd_args[1:])  # Skip script name
    else:
        # Run mm_run.main and capture returned args in dry run mode
        result = mm_run.main(cmd_args[1:])  # Passing the argument list without script name

        # Validate the returned values against expected values
        assert result.steps == expected_steps, f"Expected steps: {expected_steps}, but got: {result.steps}"
        assert result.numFrames == expected_frames, f"Expected frames: {expected_frames}, but got: {result.numFrames}"
        assert result.reportInterval == expected_interval, f"Expected interval: {expected_interval}, but got: {result.reportInterval}"

if __name__ == "__main__":
    test_numFrames_reportFrequency_handling()
