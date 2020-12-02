#-- my_code_hw01.py
#-- hw01 GEO1015.2020
#-- Pratyush Kumar
#-- 5359252
#-- Simon Pena Pereira
#-- 5391210


#-- import outside the standard Python library are not allowed, just those:
#%%
import math
import numpy as np
import scipy.spatial
import startin 
#-----

#%%
def raster_frame_creator(np_list ,cellsize):
    """returns raster as a 1d array (list of all coordinates of pixels referring to lower left)
    the no of cells as a tuple and the bounding box         

    Args:
        np_list ([numpy array]): [numpy array of list of lists sent to function]
        cellsize ([float]): [cellsize as processed and passed on to main function]
    """
    #compute bbox 
    x_list = np_list[:,0].copy()
    y_list = np_list[:,1].copy()
    z_list = np_list[:,2].copy()
    
    xmin, xmax, ymin, ymax = x_list.min(), x_list.max(), y_list.min(), y_list.max()

    #determine no of cells and create bbox
    no_x, no_y = 0, 0

    if (xmax-xmin)%cellsize == 0:
        no_x = int((xmax-xmin)//cellsize)
    else:
        no_x = int((xmax-xmin)//cellsize) +1

    if (ymax-ymin)%cellsize == 0:
        no_y = int((ymax-ymin)//cellsize)
    else:
        no_y = int((ymax-ymin)//cellsize) +1
    
    bbox = ((xmin,ymin) , (xmin + no_x*cellsize , ymin + no_y*cellsize))

    #raster creation
    rast_x = np.arange(bbox[0][0],bbox[1][0], cellsize)
    rast_y = np.arange(bbox[0][1],bbox[1][1], cellsize)
    rast_coord = np.array([[i,j] for j in rast_y for i in rast_x])
    rast_coord = np.flipud(rast_coord)
    return(rast_coord , z_list , (no_x, no_y) , bbox)

def asc_file(no_y, no_x, xmin, ymin, cellsize, filename, rast_z):
    ##writing asc file
    rast_z = rast_z.reshape(no_y , no_x)
    rast_z = np.fliplr(rast_z)
    fh = open(filename, "w")
    fh.write(f"NCOLS {no_x}\nNROWS {no_y}\nXLLCORNER {xmin}\nYLLCORNER {ymin}\nCELLSIZE {cellsize}\nNODATA_VALUE {-9999}\n") 
    for i in rast_z:
        fh.write(" ".join([str(_) for _ in i]) + '\n')
    fh.close()
    print("File written to", filename)    

#%%
def nn_interpolation(list_pts_3d, j_nn):
    """
    !!! TO BE COMPLETED !!!
     
    Function that writes the output raster with nearest neighbour interpolation
     
    Input:
        list_pts_3d: the list of the input points (in 3D)
        j_nn:        the parameters of the input for "nn"
    Output:
        returns the value of the area
 
    """  
    print("cellsize:", j_nn['cellsize'])
    cellsize =  float(j_nn['cellsize'])
    #compute bbox     #convert list3d to numpy array find min and max coordinates
    np_list = np.array(list_pts_3d)
    #make raster frame as 1d
    rast_coord , z_list,(no_x,no_y) ,bbox= raster_frame_creator(np_list , cellsize)
    xmin , ymin = bbox[0]

    list_pts = np_list[:,[0,1]]
    kd = scipy.spatial.KDTree(list_pts)
    
    rast_z = []
    dt = scipy.spatial.Delaunay(np_list[:,[0,1]])
    for coord in rast_coord:
        tri_indx =  dt.find_simplex(coord)
        if (not tri_indx) :
            # NODATA
            rast_z.append(-9999)
            continue
        elif (tri_indx == -1):
            # NODATA
            rast_z.append(-9999)            
            continue
        _ , indx = kd.query(coord, k=1)
        rast_z.append(z_list[indx])
    
    #to put in the interpolation for z values
    rast_z=np.array(rast_z)
    rast_z=rast_z

    filename = j_nn['output-file']
    asc_file(no_y, no_x, xmin, ymin, cellsize, filename, rast_z)


#%%
def idw_interpolation(list_pts_3d, j_idw):
    """
    !!! TO BE COMPLETED !!!
     
    Function that writes the output raster with IDW
     
    Input:
        list_pts_3d: the list of the input points (in 3D)
        j_idw:       the parameters of the input for "idw"
    Output:
        returns the value of the area
 
    """  
    print("cellsize:", j_idw['cellsize'])
    print("radius:", j_idw['radius'])

    cellsize =  float(j_idw['cellsize'])
    radius =  float(j_idw['radius'])
    power =  float(j_idw['power'])
    np_list = np.array(list_pts_3d)
    
    rast_coord , z_list , (no_x, no_y), bbox = raster_frame_creator(np_list ,cellsize)
    x_list = np_list[:,0].copy()
    y_list = np_list[:,1].copy()
    z_list = np_list[:,2].copy()
    xmin , ymin = bbox[0]
    
    list_pts = np_list[:,[0,1]]
    
    kd = scipy.spatial.KDTree(list_pts)
    idw_rast_z = []
    dt = scipy.spatial.Delaunay(np_list[:,[0,1]])
    for coord in rast_coord:
        tri_indx =  dt.find_simplex(coord)
        if (not tri_indx) :
            # NODATA
            idw_rast_z.append(-9999)
            continue
        elif (tri_indx == -1):
            # NODATA
            idw_rast_z.append(-9999)            
            continue

        i = kd.query_ball_point(coord, radius)
        if not i: 
            idw_rast_z.append(-9999)
        else:         
            weights = []
            known_z = []

            for indx in i:
                i_x, i_y = coord[0], coord[1] 
                p_x, p_y = list_pts[indx][0], list_pts[indx][1]
                
                if np.all(list_pts[indx] == coord):
                    weight = 1
                    z = z_list[indx]

                    weights.append(weight)
                    known_z.append(z) 
                else: 
                    dist = ((p_x - i_x)**2 + (p_y - i_y)**2)
                    weight = (1/(dist)**power)
                    z = z_list[indx]
                    
                    weights.append(weight)
                    known_z.append(z) 

            w_array = np.array(weights)
            z_array = np.array(known_z)

            z_value = (sum(w_array * z_array)/sum(w_array))
            idw_rast_z.append(z_value)

    rast_z = np.array(idw_rast_z)
    rast_z = rast_z

    filename = j_idw['output-file']

    asc_file(no_y, no_x, xmin, ymin, cellsize, filename, rast_z)
    

#%%
def tin_interpolation(list_pts_3d, j_tin):
    """
    !!! TO BE COMPLETED !!!
     
    Function that writes the output raster with linear in TIN interpolation
     
    Input:
        list_pts_3d: the list of the input points (in 3D)
        j_tin:       the parameters of the input for "tin"
    Output:
        returns the value of the area
 
    """  
    #take params and create a raster outline as in nn
    cellsize =  float(j_tin['cellsize'])
    np_list = np.array(list_pts_3d)
    # load the rast coord and the other things using predefined function
    rast_coord , z_list,(no_x,no_y) ,bbox= raster_frame_creator(np_list , cellsize)
    xmin , ymin = bbox[0]
    # rast_coord = rast_coord.reshape(int(no_x),int(no_y))
    rast_z = []
    # delauney triangulation of the x and y values obtained from file
    dt = scipy.spatial.Delaunay(np_list[:,[0,1]])

    # find triangles
    for coord in rast_coord:
        tri_indx =  dt.find_simplex(coord)
        if (not tri_indx) :
            # NODATA
            rast_z.append(-9999)
            continue
        elif (tri_indx == -1):
            # NODATA
            rast_z.append(-9999)            
            continue

        vert =  np_list[ dt.simplices[tri_indx] ,:] # coordinates of the vertices of the triangle
        #calculate barycentric weights
       
        denom = (((vert[1][1] - vert[2][1])*(vert[0][0]-vert[2][0])) + ((vert[2][0]-vert[1][0])*(vert[0][1]-vert[2][1])))

        w1_nom = ((vert[1][1] - vert[2][1])*(coord[0] - vert[2][0])) + ((vert[2][0]-vert[1][0])*(coord[1]-vert[2][1])) 
        w2_nom = ((vert[2][1] - vert[0][1])*(coord[0] - vert[2][0])) + ((vert[0][0]-vert[2][0])*(coord[1] - vert[2][1])) 
        w1= w1_nom/denom
        w2= w2_nom/denom
        w3 = 1 - w1 - w2
        # # once weight found multiply weight with the z values at vertex of each triangle
        z_val = (vert[0][2]*w1) + (vert[1][2]*w2) + (vert[2][2]*w3)
        # print(z_val)
        rast_z.append(z_val)
    #write to file
    rast_z = np.array(rast_z)#.reshape(no_x, no_y)    
    filename = j_tin['output-file']
    asc_file(no_y, no_x, xmin, ymin, cellsize, filename, rast_z)
    
#%%
def kriging_interpolation(list_pts_3d, j_kriging):
    """
    !!! TO BE COMPLETED !!!
     
    Function that writes the output raster with ordinary kriging interpolation
     
    Input:
        list_pts_3d: the list of the input points (in 3D)
        j_kriging:       the parameters of the input for "kriging"
    Output:
        returns the value of the area
 
    """  
    
    cellsize =  float(j_kriging['cellsize'])
    radius =  float(j_kriging['radius'])

    np_list = np.array(list_pts_3d)
    
    rast_coord , z_list , (no_x, no_y), bbox = raster_frame_creator(np_list ,cellsize)

    z_list = np_list[:,2].copy()
    xmin , ymin = bbox[0]
    
    list_pts = np_list[:,[0,1]]
    dt_start = startin.DT()
    dt_start.insert(list_pts_3d)
    dt_no_dup = dt_start.all_vertices()[1:]
    dt_pt_2d = [(x,y) for x,y,z in dt_no_dup]
    z_list = [(z) for x,y,z in dt_no_dup]
    kd = scipy.spatial.KDTree(dt_pt_2d, 50)
    krig_rast_z = []
    # from the variogram.py file
    # using gaussian
    nugget = 1
    sill = 1375
    rng = 285

    gam = lambda y : (nugget +sill * (1.0 - math.exp(-9.0*y*y/(rng**2))))
    # with i j as points of type x,y
    dist_xy = lambda i,j : (math.sqrt((j[1] - i[1])**2 + (j[0] - i[0])**2) )
    cnt=0

    dt = scipy.spatial.Delaunay(np_list[:,[0,1]])
        # create dt

    for coord in rast_coord:

        tri_indx =  dt.find_simplex(coord)
        if (not tri_indx) :
            # NODATA
            krig_rast_z.append(-9999)
            continue
        elif (tri_indx == -1):
            # NODATA
            krig_rast_z.append(-9999)            
            continue
        if not dt_start.locate(*coord):
            krig_rast_z.append(-9999)
            continue
        i = kd.query_ball_point(coord, radius)
        if not i: 
            krig_rast_z.append(-9999)
            continue

        neighbor_coords = [dt_pt_2d[indx] for indx in i]
        z_neighbor = [z_list[indx] for indx in i]
        dist_from_neighbors = [dist_xy(neighbor , coord) for neighbor in neighbor_coords]
        
        #when one of the points is a neighbor of itself, or duplicate point
        if 0 in dist_from_neighbors:
            i = dist_from_neighbors.index(0)
            krig_rast_z.append(z_neighbor[i])
            continue
        
        #creating a covariance matrix for variogram
        cov = []
        for i in neighbor_coords:
            cov_row = []
            for j in neighbor_coords:
                cov_row.append(gam(dist_xy(i,j)))
            cov.append(cov_row+[1])
            # cov.append([1]) #covariance matrix has 1 at end
        cov.append([1]*len(neighbor_coords)+[0]) # last row
        cov = np.array(cov) #convert to numpy array
        
        #creatnig d matrix
        d=[]
        for i in neighbor_coords:
            d.append([gam(dist_xy(i, coord))])
        d.append([1]) # last value in vector
        d=np.array(d ) #convert to np array

        try:
            #cov inverse @ d
            w = np.matmul(np.linalg.inv(cov) , d)
        except:
            krig_rast_z.append(-9999)

            continue
        #calculate weight matrix
        weights_arr = [weight_val[0] for weight_val in w[:-1]]
        
        normalized_w = [w_/sum(weights_arr) for w_ in weights_arr] if sum(weights_arr) else weights_arr
        z = sum([z_*wt for z_ ,wt in zip(z_neighbor, normalized_w)])
        krig_rast_z.append(z)

    rast_z = np.array(krig_rast_z)
    rast_z = rast_z

    filename = j_kriging['output-file']
    asc_file(no_y, no_x, xmin, ymin, cellsize, filename, rast_z)
# %%
