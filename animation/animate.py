import multiprocessing
import subprocess as sp
import tqdm
import os

#TODO Add info to file concerning resolution

#   Settings

F_START = 0
F_END = 1
FRAMES = 600
FRAME_RATE = 60
MAGNETS = 3



FFMPEG_BIN = "ffmpeg"

CMD_COMPILE = [
    "gcc",
    "fast_magnet.c",
    "-o",
    "magnet",
    "-lm"    
]

CMD_MAGNET = [
    "magnet",
    "-q",
    "-f", "%lg",
    "-m", "%d", 
    "-o", "%s"
]

CMD_FFMPEG = [
    FFMPEG_BIN,
    "-loglevel", "info",
    "-s", "1920x1080",
    "-i", "_frames//img%04d.png",
    "test.avi"
]

def _foo(i):
    pngFile = "_frames//img%04d.png" % i
    force = (F_END - F_START) / FRAMES * i

    cmd = " ".join(CMD_MAGNET) % (force, MAGNETS, pngFile)

    if not os.path.isfile(pngFile):
        process = sp.Popen(cmd, shell=True, stdout=sp.PIPE)
        process.wait()
        if process.returncode != 0:
            raise(Exception)



if __name__ == '__main__':

    sp.run(" ".join(CMD_COMPILE))
    print("successfully compiled code")

    with multiprocessing.Pool(multiprocessing.cpu_count()) as p:
        r = list(tqdm.tqdm(p.imap_unordered(_foo, range(1, FRAMES + 1)), total=FRAMES))

    sp.run(" ".join(CMD_FFMPEG))

