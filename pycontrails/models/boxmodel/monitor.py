import psutil
import os
import time

# Get process info
pid = os.getpid()
py = psutil.Process(pid)

while True:
    # Print CPU info
    cpu_usage = os.popen("ps aux | grep " + str(pid)).read().split()[2]
    print("CPU usage : {}%".format(cpu_usage))

    # Print memory info
    memory_usage = py.memory_info()[0]/2.**30  # memory use in GB
    print("Memory usage : {} GB".format(memory_usage))
    
    # Sleep for a while
    time.sleep(0.5)