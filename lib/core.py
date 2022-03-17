#############################################################################
# Core functions for the degenotate program
# Gregg Thomas
#############################################################################

import sys
import os
import timeit
import datetime
import subprocess

#############################################################################

def errorOut(errnum, errmsg, globs):
# Formatting for error messages.
    fullmsg = "**Error " + str(errnum) + ": " + errmsg;
    border = "-" * len(fullmsg);
    fullstr = "\n" + border + "\n" + fullmsg + "\n" + border + "\n"
    printWrite(globs['logfilename'], globs['log-v'], "\n" + border + "\n" + fullmsg + "\n" + border + "\n");
    if globs['endprog']:
        globs['exit-code'] = 1;
        endProg(globs);
    else:
        printWrite(globs['logfilename'], globs['log-v'], "\nScript call: " + " ".join(sys.argv));
        sys.exit(1);

#############################################################################

def fileCheck(globs):
# Checks file options.
    files = ['gxf-file', 'fa-file', 'in-seq'];
    for f in files:
        if globs[f]:
            if not os.path.isfile(globs[f]) and not os.path.isdir(globs[f]):
                errorOut("CORE1", "File/path not found: " + globs[f], globs);
            globs[f] = os.path.abspath(globs[f]);
    return globs;

#############################################################################

# def execCheck(globs, a):
# # Checks dependency executables.
#     deps_passed = True;
#     # Variable to check if all dependencies are found.

#     if a.phyloacc_path:
#         globs['phyloacc'] = a.phyloacc_path;
#     # Update the global paths if the user provided them through args.

#     dpad = 20;
#     if a.depcheck:
#         print("# --depcheck set: CHECKING DEPENDENCY PATHS AND EXITING.\n");
#         print(spacedOut("   PROGRAM", dpad) + spacedOut("PATH", dpad) + "STATUS");
#         print("   -------------------------------------------");
#     # For the dependency check option (--depcheck), this initializes a neat output table.

#     for opt in ['phyloacc', 'phyloacc-gbgc']:
#         dcheck_str = [spacedOut("   " + opt, dpad), spacedOut(globs[opt], dpad), "NA"];
#         # Initialize the check string for --depcheck.

#         cmd_result = subprocess.run(globs[opt], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE);
#         # Run the provided command and retrieve the exit code.

#         if cmd_result.returncode > 1:
#         # If the exit code for the command run is greater than 1, the command isn't found.
#             dcheck_str[2] = "FAILED with exit code " + str(cmd_result.returncode);
#             deps_passed = False;
#             # Update the check string and keep going.
#             if not a.depcheck:    
#                 errorOut("CORE2", prog + " not found at specified path: " + globs[opt], globs);
#             # On a normal run, exit immediately.
#         else:
#             dcheck_str[2] = "PASSED";
#             # Update the check string.
            
#         if a.depcheck:
#             print("".join(dcheck_str));
#         # Print the check string if --depcheck is set.

#     return globs, deps_passed;

#############################################################################

def detectCompression(filename):
# Detect compression of a file by examining the first lines in the file

    compression_type = "none";

    magic_dict = {
            b"\x1f\x8b\x08": "gz",
            # b"\x1f\x8b\x08\x08": "gz",
            b"\x42\x5a\x68": "bz2",
            b"\x50\x4b\x03\x04": "zip"
        }
    # An encoded set of possible "magic strings" that start different types of compressed files
    # From: https://www.garykessler.net/library/file_sigs.html
    # \x is the escape code for hex values
    # b converts strings to bytes

    max_len = max(len(x) for x in magic_dict)
    # The number of characters to read from the beginning of the file should be the length of
    # the longest magic string

    file_start = open(filename, "rb").read(max_len);
    # Read the beginning of the file up to the length of the longest magic string

    for magic_string in magic_dict:
        if file_start.startswith(magic_string):
            compression_type = magic_dict[magic_string];
    # Check each magic string against the start of the file

    return compression_type;

#############################################################################

def getOutTime():
# Function to get the date and time in a certain format.
    return datetime.datetime.now().strftime("%m-%d-%Y.%I-%M-%S");

#############################################################################

def getDate():
# Function to get the date and time in a certain format.
    return datetime.datetime.now().strftime("%m.%d.%Y");

#############################################################################

def getTime():
# Function to get the date and time in a certain format.
    return datetime.datetime.now().strftime("%H:%M:%S");

#############################################################################

def getDateTime():
# Function to get the date and time in a certain format.
    return datetime.datetime.now().strftime("%m.%d.%Y | %H:%M:%S");

#############################################################################

def isPosInt(numstr):
# Check if a string is a positive integer
    try:
        num = int(numstr);
    except:
        return False;

    if num > 0:
        return num;
    else:
        return False;

#############################################################################

def printWrite(o_name, v, o_line1, o_line2="", pad=0):
# Function to print a string AND write it to the file.
    if o_line2 == "":
        outline = o_line1;
    else:
        outline = o_line1 + " "*(pad-len(o_line1)) + o_line2;
    if v in [-1,1,2]:
        print(outline);
    if v != -1:
        f = open(o_name, "a");
        f.write(outline + "\n");
        f.close();

#############################################################################
    
def spacedOut(string, totlen, sep=" "):
# Properly adds spaces to the end of a message to make it a given length
    spaces = sep * (totlen - len(string));
    return string + spaces;

#############################################################################

def report_step(globs, step, step_start_time, step_status, start=False, full_update=False):
# Uses psutil to gather memory and time info between steps and print them to the screen.

    dashes = 150
    if globs['psutil']:
        import psutil;
        dashes = 175;
    # Determine the number of dashes to frame the update table depending on the presence of psutil

    cur_time = timeit.default_timer();
    # The time at the start of the status update

    col_widths = [ 14, 10, 50, 40, 20, 16 ];
    if globs['psutil']:
        col_widths += [18, 10];
    # The column widths

    if start:
        headers = [ "# Date", "Time", "Current step", "Status", "Elapsed time (s)", "Step time (s)" ];
        if globs['psutil']:
            headers += ["Current mem (MB)", "Virtual mem (MB)"]
        # A list of the headers

        headers = "".join([ spacedOut(str(headers[i]), col_widths[i]) for i in range(len(headers)) ]);
        # Converting the list to a string based on the column widths

        printWrite(globs['logfilename'], globs['log-v'], "# " + "-" * dashes);
        printWrite(globs['logfilename'], globs['log-v'], headers);
        printWrite(globs['logfilename'], globs['log-v'], "# " + "-" * dashes);
        # Print the dashes and the headers
    # The first call is just to print the headers

    ##########

    else:
        prog_elapsed = str(round(cur_time - globs['starttime'], 5));
        # Get the total amount of time that the program has been running

        if not step_start_time:
        # If no step start time is given, then this is the first entry for this status
        # update, that will display "In progress..." or similar.

            out_line = [ "# " + getDate(), getTime(), step, step_status ];
            # The output for the initial status entry includes the date, time, step label, and progress message

            term_col_widths = col_widths[:4];
            # Only get the first 4 column widths for the initial status entry

            out_line = [ spacedOut(str(out_line[i]), term_col_widths[i]) for i in range(len(out_line)) ];

            if full_update:
                out_line += "\n";
            # For some status updates, intermediate info will be printed, in which case we add a newline here

            if not globs['quiet']:
                sys.stdout.write("".join(out_line));
                sys.stdout.flush();
            # Convert the output list to a string, write, and flush stdout

        ## The initial status entry to display "In progress..."
        #####

        else:
            step_elapsed = str(round(cur_time - step_start_time, 5));
            # Get the full step time here

            out_line = [ step_status, prog_elapsed, step_elapsed ];
            # Gather info for the full output line to print to screen

            if globs['psutil']:
                mem = round(sum([p.memory_info()[0] for p in globs['pids']]) / float(2 ** 20), 5);
                vmem = round(sum([p.memory_info()[1] for p in globs['pids']]) / float(2 ** 20), 5);
                out_line += [str(mem), str(vmem)];
            # If psutil is present, get current memory info

            term_col_widths = col_widths[3:];
            # Get the column widths for the print to screen output

            file_line = [ "# " + getDate(), getTime(), step ] + out_line;
            file_col_widths = col_widths[:3] + [30] + col_widths[4:];
            # For output to the file, we write the whole line each time
            # Add the initial entry fields here
            # This will also be used for some status updates where the whole message needs to be printed
            # to the screen
            
            out_line = [ spacedOut(str(out_line[i]), term_col_widths[i]) for i in range(len(out_line)) ];
            file_line = [ spacedOut(str(file_line[i]), col_widths[i]) for i in range(len(file_line)) ];
            # Compile both the truncated and the full status update

            if not globs['quiet']:
                if full_update:
                    sys.stdout.write("".join(file_line) + "\n");
                    sys.stdout.flush();
                else:         
                    sys.stdout.write("\b" * 40);
                    sys.stdout.write("".join(out_line) + "\n");
                    sys.stdout.flush();
            # For full updates, print the full line to the screen
            # For others, delete the "In progress..." column and update the same status line
            
            printWrite(globs['logfilename'], 0, "".join(file_line));
            # Write the full line to the file.
        # The final status entry
        #####

    return cur_time;

#############################################################################

def welcome():
# Reads the ASCII art "Referee" text to be printed to the command line.
    return open(os.path.join(os.path.dirname(__file__), "welcome.txt"), "r").read();

#############################################################################

def endProg(globs):
# A nice way to end the program.
    if globs['quiet']:
        globs['log-v'] = 1;
    endtime = timeit.default_timer();
    totaltime = endtime - globs['starttime'];

    printWrite(globs['logfilename'], globs['log-v'], "# " + "=" * 175);
    printWrite(globs['logfilename'], globs['log-v'], "#\n# Done!");
    printWrite(globs['logfilename'], globs['log-v'], "# The date and time at the end is: " + getDateTime());
    printWrite(globs['logfilename'], globs['log-v'], "# Total execution time:            " + str(round(totaltime,3)) + " seconds.");
    printWrite(globs['logfilename'], globs['log-v'], "# Output directory for this run:   " + globs['outdir']);
    printWrite(globs['logfilename'], globs['log-v'], "# Log file for this run:           " + globs['logfilename']);

    # if globs['aln-stats-written']:
    #     printWrite(globs['logfilename'], globs['log-v'], "# Alignment stats file:            " + globs['alnstatsfile']);    

    # if globs['scf-stats-written']:
    #     printWrite(globs['logfilename'], globs['log-v'], "# Concordance factor stats file:   " + globs['scfstatsfile']); 

    # if globs['scf-tree-written']:
    #     printWrite(globs['logfilename'], globs['log-v'], "# Concordance factor tree file:    " + globs['scftreefile']); 

    if globs['exit-code'] != 0:
        printWrite(globs['logfilename'], globs['log-v'], "#\n# ERROR: NON-ZERO EXIT STATUS.");
        printWrite(globs['logfilename'], globs['log-v'], "# ERROR: DEGENOTATE FINISHED WITH ERRORS.");
        printWrite(globs['logfilename'], globs['log-v'], "# ERROR: PLEASE CHECK THE LOG FILE FOR MORE INFO: " + globs['logfilename'] + "\n#");

    if globs['warnings'] != 0:
        printWrite(globs['logfilename'], globs['log-v'], "\n# degenotate finished with WARNINGS -- check log file for more info");

    #print("# " + "=" * 125);
    printWrite(globs['logfilename'], globs['log-v'], "# " + "=" * 175);
    printWrite(globs['logfilename'], globs['log-v'], "#");
    sys.exit(globs['exit-code']);

#############################################################################