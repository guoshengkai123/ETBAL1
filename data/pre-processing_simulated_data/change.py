def Get_file_Camera_parameter_data(file_name):

    parameter_list = []
    parameter_file = open(file_name, 'r')
    for line in parameter_file:
        if 'C' not in line:
            line = line.strip('\n')
            line = line.split('\t')
            parameter_list.append(line)
    for item in parameter_list:
        for i in range(0, len(item)):
            item[i] = float(item[i])

    a=[]
    for i in range(len(parameter_list)):
        for j in range(6):
            a.append(parameter_list[i][j])

    return a



def Get_file_Points_data(file_name):

    points_list = []
    points_file = open(file_name, 'r')
    for line in points_file:
        line = line.strip('\n')
        line = line.strip('[')
        line = line.strip(']')
        line = line.split(', ')
        points_list.append(line)
    for item in points_list:
        for i in range(0, len(item)):
            item[i] = float(item[i])
    a=[]
    for i in range(len(points_list)):
        for j in range(3):
            a.append(points_list[i][j])


    return a


def Get_file_Points3d_data(file_name,ncams):

    points_list = []
    points_file = open(file_name, 'r')
    for line in points_file:
        line = line.replace(' 0,', ' 0.0, 0.0,')
        line = line.replace('[0,', '[0.0, 0.0,')
        line = line.replace(' 0]', ' 0.0, 0.0]')
        line = line.replace('[', '')
        line = line.replace(']', '')
        line = line.replace(',', ' ')
        line = line.strip('\n')
        line = line.split()
        points_list.append(line)

    for item in points_list:
        for i in range(0, len(item)):
            item[i] = float(item[i])
    a=[]
    for i in range(len(points_list)):
        for j in range(ncams*2):
            a.append(points_list[i][j])

    #print(points_list)
    return a


def Store_txt(data, file_name):

    file = open(file_name, 'w')
    for i in range(len(data)):
        file.write(str(data[i]))
        file.write('\n')
    file.close()

def Store_mask(length, file_name):

    file = open(file_name, 'w')
    for i in range(length):
        file.write('1')
        file.write('\n')
    file.close()

if __name__ == '__main__':
    cam_list=Get_file_Camera_parameter_data('64_cameras.txt')
    #print(cam_list)
    Store_txt(cam_list, 'changedata_output/cams.txt')

    pots3D=Get_file_Points_data('40_points.txt')
    #print(pots3D)
    Store_txt(pots3D, 'changedata_output/points3d.txt')

    pots2D=Get_file_Points3d_data('40points_64cameras_Proj_GroundTruth_3to2.txt',(int)(len(cam_list)/6))
    #print(pots2D)
    Store_txt(pots2D, 'changedata_output/points2d.txt')

    Store_mask((int)(len(pots2D)/2),'changedata_output/vmask.txt')
    
