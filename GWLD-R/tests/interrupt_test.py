"""Check that R can cancel native matrix and streamed pair computations."""
import os
import pathlib
import shutil
import signal
import subprocess
import tempfile
import time


def check(mode):
    with tempfile.TemporaryDirectory(prefix="gwld-interrupt-") as directory:
        root = pathlib.Path(directory)
        script = root / "interrupt.R"
        script.write_text('''library(GWLD)
args <- commandArgs(trailingOnly=TRUE)
set.seed(2026)
g <- matrix(sample(0:2, 256L*6000L, replace=TRUE), 256L, 6000L)
info <- cbind(rep("1",6000L),seq_len(6000L),paste0("v",seq_len(6000L)))
writeLines("ready",args[2])
tryCatch({
  if (args[1]=="matrix") GWLD(g,method="RMI")
  else calc_decay(g,info,method="RMI",maxdist=300L,file=args[3])
}, interrupt=function(e) quit(status=42L,save="no"))
quit(status=1L,save="no")
''')
        marker = root / "ready"
        output = root / "pairs.tsv"
        process = subprocess.Popen(
            [shutil.which("Rscript"), "--vanilla", str(script), mode, str(marker), str(output)],
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
        )
        try:
            deadline = time.monotonic() + 30
            while not marker.exists():
                if process.poll() is not None or time.monotonic() > deadline:
                    raise AssertionError("R computation did not start: " + process.communicate()[0])
                time.sleep(0.02)
            time.sleep(0.1)
            os.kill(process.pid, signal.SIGINT)
            text, _ = process.communicate(timeout=10)
            if process.returncode != 42:
                raise AssertionError(f"{mode} did not handle interrupt: {process.returncode}\n{text}")
            print(f"{mode}: native computation interrupted successfully")
        finally:
            if process.poll() is None:
                process.kill()
                process.communicate()


if __name__ == "__main__":
    check("matrix")
    check("stream")
