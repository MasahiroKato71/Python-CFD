from __future__ import annotations


def MaxAbs(*args):
    maxNum, minNum = max(args), min(args)
    return maxNum if abs(maxNum) >= abs(minNum) else minNum