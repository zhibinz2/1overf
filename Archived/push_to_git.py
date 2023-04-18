# Read and write files using the built-in Python file methods
import subprocess

def main():  
    run("add","-A")
    run("commit","-am","update")
    # run("push", "-u", "origin", "master")   

def run(*args):
    return subprocess.check_call(['git'] + list(args))
    
if __name__ == "__main__":
    main()