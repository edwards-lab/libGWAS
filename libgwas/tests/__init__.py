
import os
import libgwas

def remove_file(fn):
    try:
        libgwas.close_file(fn)
        os.remove(fn)
    except:
        sys.stderr.write(f"Unable to close file, '{fn}'")
        
