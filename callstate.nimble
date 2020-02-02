# Package

version       = "0.0.2"
author        = "Luobin Yang"
description   = "GATK3 CallableLoci Replacement"
license       = "MIT"

# Dependencies

requires "hts >= 0.3.2", "nim >= 1.0.0", "docopt >= 0.6.8"

srcDir = "src"

bin = @["callstate"]
skipDirs = @["tests"]


task test, "run the tests":
  exec "nim c --lineDir:on --debuginfo -r tests/test_callstate.nim"
