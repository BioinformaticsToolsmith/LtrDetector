import os

def buildDirectories():

    os.mkdir("../bin")
    os.mkdir("../bin/exception")
    os.mkdir("../bin/nonltr")
    os.mkdir("../bin/test")
    os.mkdir("../bin/tr")
    os.mkdir("../bin/utility")
    
if __name__ == "__main__":

    buildDirectories()