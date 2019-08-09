
import os

def remove_file(fn):
    try:
        close_file(fn)
        os.remove(fn)
    except:
        sys.stderr.write(f"Unable to close file, '{fn}'")
        
def close_file(f):
    '''Windows...I imagine this requires some other changes, but until
    the rest of this is working, this will eliminate the stupid failures'''
    
    if f is not None:
        try:
            f.close()
        except:
            pass
