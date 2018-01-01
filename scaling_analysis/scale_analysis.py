import numpy as np
import matplotlib.pyplot as plt
import matplotlib

def scale_analysis(cpu_file, gpu_file, ns, title):
    cpus = np.loadtxt(cpu_file)
    gpus = np.loadtxt(gpu_file)
    plt.loglog(ns, cpus, 'b-o', basex=2, basey=2, label="cpu time")
    plt.loglog(ns, gpus, 'r--o', basex=2, basey=2, label="gpu time")
    plt.legend(loc="upper left");
    plt.xlabel('number of particles')
    plt.ylabel('times (s)')
    plt.title(title)
    plt.show()

if __name__ == "__main__":
    matplotlib.rcParams.update({'font.size': 22})
    '''
    cpu_file = "dynamic_times_cpu.txt"
    gpu_file = "dynamic_times_gpu.txt"
    title = "Dynamic Simulation"
    ns = np.array([2, 4, 10, 20, 50, 100])
    '''
    cpu_file = "thermal_times_cpu.txt"
    gpu_file = "thermal_times_gpu.txt"
    title = "Thermal Simulation"
    lst = [2**x for x in range(7, 15)]
    ns = np.array(lst)
    scale_analysis(cpu_file, gpu_file, ns, title)
