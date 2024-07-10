import numpy as np


def write_evpfft_input(fname, grain_ids, orientation, nx_ny_nz):
    count = 0
    f = open(fname, 'w')
    for z in range(nx_ny_nz[2]):
        for y in range(nx_ny_nz[1]):
            for x in range(nx_ny_nz[0]):
                my_line=''
                for i in range(3):
                    my_line += str(orientation[grain_ids[count]-1,i])
                    my_line += " "
                my_line += str(x+1) + " " + str(y+1) + " " + str(z+1) + " " + str(grain_ids[count]) + " " + str(1) + "\n"
                f.write(my_line)
                count+=1
    f.close()
    return
def read_neper_raster(fname):
    file1 = open(fname, 'r')
    while True:

        # Get next line from file
        line = file1.readline().strip()

        # if line is empty
        # end of file is reached
        if not line:
            break
        if line == "**general":
            dims = int(file1.readline())
            nx_ny_nz = np.asarray(file1.readline().strip().split())
            nx_ny_nz = nx_ny_nz.astype(int)
            nvoxel = np.prod(nx_ny_nz)
        if line == "**cell":
            ngrains = int(file1.readline().strip())
            ori = np.zeros([ngrains,3],dtype=float)
            ngrains_count = 0
            dummy = file1.readline().strip()
            while(True) :
                line =  file1.readline().strip()
                if line[0]=="*":
                    break
                grain_ids_temp = np.asarray(line.split())
                grain_ids_temp = grain_ids_temp.astype(int)
                ngrains_count += np.size(grain_ids_temp)
                if ngrains_count == ngrains :
                    break
            # grain_ids = list(map(int, file1.readline().strip().split(" ")))
        if line == "*ori":
            dummy = file1.readline().strip()
            for i in range(ngrains):
                ori_temp = np.asarray(file1.readline().strip().split())
                ori[i,:] = ori_temp.astype(float)

        if line == "**data":
            grain_ids = np.zeros([nvoxel],dtype=int,order='F')
            n_read_gid = 0
            while (n_read_gid < nvoxel):
                gid_temp = np.asarray(file1.readline().strip().split())
                gid_temp = gid_temp.astype(int)
                n_read_temp = np.size(gid_temp)
                grain_ids[n_read_gid:(n_read_gid+n_read_temp)]=gid_temp
                n_read_gid += n_read_temp

    file1.close()
    return grain_ids, ori, nx_ny_nz

if __name__ == "__main__":
    import sys
    fname_tesr = sys.argv[1]+".tesr"
    fname_fft = sys.argv[1]+".dat"
    # fname = "16Cube10Grains.tesr"
    # fname_fft = "16Cube10Grains.dat"
    print("fname", fname_tesr)
    print("fname_fft", fname_fft)
    grain_ids, orientation, nx_ny_nz = read_neper_raster(fname_tesr)
    write_evpfft_input(fname_fft, grain_ids, orientation, nx_ny_nz)
