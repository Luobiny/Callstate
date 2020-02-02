import ../callstate
import unittest

suite "test procedures in callstate":
  setup:
    var raw: seq[int32] =  @[0i32, 30i32, 30i32, 16i32, 8i32, 30i32, 30i32, 40i32, 5000i32, 5i32]

  test "mean_depth":
    assert mean_depth(raw[0 ..< 2]) == 15.0
    assert mean_depth(raw[0 ..< 3]) == 20.0
    assert mean_depth(raw[0 .. 4]) == mean_depth(raw[0 ..< 5])
    assert mean_depth(raw[0 .. 7]) == 23.0 
