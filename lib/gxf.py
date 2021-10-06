#############################################################################
# Functions to handle gtf and gff annotation files.
#############################################################################

import sys
import os
import gzip
import lib.core as CORE

#############################################################################

def read(globs):

    step = "Detecting compression of annotation file";
    step_start_time = CORE.report_step(globs, step, False, "In progress...");
    globs['gxf-compression'] = CORE.detectCompression(globs['gxf-file']);
    if globs['gxf-compression'] == "none":
        reader = open;
        readline = lambda l : l.strip().split("\t");
        step_start_time = CORE.report_step(globs, step, step_start_time, "Success: No compression detected");
    else:
        reader = gzip.open;
        readline = lambda l : l.decode().strip().split("\t");
        step_start_time = CORE.report_step(globs, step, step_start_time, "Success: " + globs['gxf-compression'] + " detected");
    # Detect the compression of the input sequence file

    globs['annotation'] = {};


    # step = "Reading genes";
    # step_start_time = CORE.report_step(globs, step, False, "In progress...");

    # for line in reader(globs['gxf-file']):
    #     line = readline(line);


    return globs;