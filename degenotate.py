#!/usr/bin/env python3
#############################################################################
# Degenotate is a script to calculate degeneracy of coding sites within a 
# genome. This is the main interface.
#
# Gregg Thomas
# Summer 2021
#############################################################################

import sys
import os
import lib.core as CORE
import lib.params as params
import lib.opt_parse as OP
import lib.gxf as gxf
# import lib.seq as SEQ
# import lib.tree as TREE
# import lib.output as OUT
# import lib.batch as BATCH

#############################################################################

if __name__ == '__main__':
# Main is necessary for multiprocessing to work on Windows.

    globs = params.init();
    # Get the global params as a dictionary.
    
    print("\n" + " ".join(sys.argv) + "\n");

    if any(v in sys.argv for v in ["--version", "-version", "--v", "-v"]):
        print("# PhyloAcc version " + globs['interface-version'] + " released on " + globs['releasedate'])
        sys.exit(0);
    # The version option to simply print the version and exit.
    # Need to get actual PhyloAcc version for this, and not just the interface version.

    print("#");
    print("# " + "=" * 125);
    print(CORE.welcome());
    if "-h" not in sys.argv:
        print("            Degeneracy annotation of transcripts\n");
    # A welcome banner.

    globs = OP.optParse(globs);
    # Getting the input parameters from optParse.

    if globs['info']:
        print("# --info SET. EXITING AFTER PRINTING PROGRAM INFO...\n#")
        sys.exit(0);
    if globs['norun']:
        print("# --norun SET. EXITING AFTER PRINTING OPTIONS INFO...\n#")
        sys.exit(0);
    # Early exit options

    step_start_time = CORE.report_step(globs, "", "", "", start=True);
    # Initialize the step headers

    # if globs['gxf-file']:
    #    globs = gxf.read(globs);

    # globs = SEQ.read(globs);
    # # Library to read input sequences


    CORE.endProg(globs);

#############################################################################

