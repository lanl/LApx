import numpy as np


def write_evpfft_input(fname, grain_ids, orientation, nx_ny_nz):
    count = 0
    f = open(fname, 'w')
    f.write("LApx-wPF\n")
    for z in range(nx_ny_nz[2]):
        for y in range(nx_ny_nz[1]):
            for x in range(nx_ny_nz[0]):
                my_line=''
                my_line += str(x+1) + " " + str(y+1) + " " + str(z+1) + " " + str(grain_ids[0,count]) + " " + str(1.)
                for i in range(3):
                    my_line += " "
                    my_line += str(orientation[i,count])
                my_line += "\n"
                f.write(my_line)
                count+=1
                my_line += "\n"
    f.close()
    return


def read_old_LApx_format(fname, dims):
    file1 = open(fname, 'r')
    counter = 0
    grain_ids = np.zeros([1,np.product(dims)],dtype=int)    
    orientation = np.zeros([3,np.product(dims)])
    phase = np.zeros([1,np.product(dims)],dtype=int)
    while counter < np.product(dims):

        # Get next line from file
        line = file1.readline().strip().split()
        grain_ids[0,counter] = line[6]
        for i in range(3):
            orientation[i,counter] = float(line[i])
        phase[0,counter] = line[7]
        counter += 1

    file1.close()
    return grain_ids, orientation, phase

if __name__ == "__main__":
    import sys
    fname_old = sys.argv[1]+".dat"
    fname_new = sys.argv[1]+"-LApx-wPF.dat"
    nx= int(sys.argv[2])
    ny= int(sys.argv[3])
    nz= int(sys.argv[4])
    dims = [nx,ny,nz]
    # fname = "16Cube10Grains.tesr"
    # fname_fft = "16Cube10Grains.dat"
    print("fname_old", fname_old)
    print("fname_new", fname_new)
    grain_ids, orientation, nx_ny_nz = read_old_LApx_format(fname_old, dims)
    write_evpfft_input(fname_new, grain_ids, orientation, dims)
