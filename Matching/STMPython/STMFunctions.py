# 2017-04-12 included [0,0,0] in the default neighbours should not be too critical, just a backup measure.
import sys
import math
import collections as col
import numpy as np
import itertools as it
import copy
import datetime
#import scipy.spatial as sps

# Take a position p and a bunch of neighbours, to create p+n for each n in neighbours
def expandneighbours(p,neighbours):
    return(list(list(map(lambda x,y: x + y, p, k)) for k in neighbours))

# Take positions p and a bunch of neighbours, to create p+n for each n in neighbours and for each p
def expandallneighbours(ps,neighbours):
    return(joinlists(list(expandneighbours(p,neighbours) for p in ps)))

# 'Flatten' a list of lists to a single list
def joinlists(lst):
    return(list([value for sublst in lst for value in sublst]))

# Custom sign function that gives either -1 or 1, also for input value 0
def Sgn(x):
    if (x<0):
        return(-1)
    else:
        return(1)

# Division that spawns -infinity in case the denominator is 0
def SpecialDivision(a,b):
    if(b==0):
        return(-math.inf)
    else:
        return(a/b)

# Mathematica style map with level specification
def maplevel(f, item, level):
    if level == 0:
        return f(item)
    else:
        return [maplevel(f, i, level - 1) for i in item]

# Return unique elements of a list, works without second argument for simple elements. Lists-like elements require a lambda function lambda x:tuple(x)
def uniqify(seq, idfun=None):
    if idfun is None:
        def idfun(x): return x
    seen = {}
    result = []
    for item in seq:
        marker = idfun(item)
        if marker in seen: continue
        seen[marker] = 1
        result.append(item)
    return result

# Gives index in which 'bin' the value n will fall with boundaries 'lst', -1 = outside
# lst has to be sorted to make this work properly
# Implemented with O(Log(N)) scaling, basic idea is to reduce by factors of 2 every time until 1 bin is left
# Greedy: it will take the first bin if the right boundary exactly matches n
def PositionSorted(boundaries,n):
    lst = list(boundaries)
    lst[0] -= 10**-8
    lst[-1] += 10**-8
    mn = 0
    mx = len(lst)-1
    if (lst[mn] <= n <= lst[mx]):
        while (mx - mn > 1):
            #print('min',mn,'max',mx)
            trial = round((mn+mx)/2)          # Banker's rounding but does not matter, still O(Log(N)) scaling
            if (n > lst[trial]):
                mn = trial
            else:
                mx = trial
        #print('Final: min',mn,'max',mx)        
        return(mn)
    else:
        return(-1)

# Euclidean distance of a list (faster than np.linalg.norm(y)!)
def VectorNorm(y):
    x = np.array(y)
    return(np.sqrt(x.dot(x)))

def SquareVectorNorm(y):
    x = np.array(y)
    return(x.dot(x))

def normalize(v):
    norm = VectorNorm(v)
    if norm==0:
        return(v)
    return(v/norm)

def ClosestPointToLines2(p,v):
    p1 = np.array(p[0])
    p2 = np.array(p[1])
    v1 = np.array(v[0])
    v2 = np.array(v[1])
    #a = np.dot(v1,v1)          # Assuming v is normalized => length 1    
    b = 2*np.dot(p1-p2,v1)
    c = 2*np.dot(v1,v2)
    d = 2*np.dot(p2-p1,v2)
    #e = np.dot(v2,v2)          # Assuming v is normalized => length 1               
    f = np.dot(p1,p1) + np.dot(p2,p2)
    #s = (2*a*d + b*c)/(c**2-4*a*e)
    s = (2*d + b*c)/(c**2-4)    # Assuming v is normalized => a=e=1
    #t = (c*s - b)/(2*a)
    t = (c*s - b)/2             # Assuming v is normalized => a=e=1
    sol = (p1 + t*v1 + p2 + s*v2)/2
    #d1 = VectorNorm(np.cross(v1,p1-sol))/np.sqrt(a)  # Assuming v is normalized => a=1
    d1 = VectorNorm(np.cross(v1,p1-sol))
    #d2 = np.linalg.norm(np.cross(sol-p1,sol-p1+v1))/np.linalg.norm(v1)  # Must be the same as d1 for two lines!
    return([sol.tolist(),d1.item()])

def ClosestPointToLines(p,v):
    if(len(p)==2):
        return(ClosestPointToLines2(p,v))
    else:
        a = np.array(p)
        #d = np.array(list(map(normalize,v)))
        d = np.array(v)   # Assuming v is normalized already
        length = len(p)
        rhs = np.array([0.0,0.0,0.0])
        lhs = length*np.identity(3)
        for i in range(length):
            rhs += a[i]-d[i]*np.dot(a[i],d[i])
            lhs -= np.outer(d[i],d[i])
        sol = np.linalg.solve(lhs, rhs)
        dists = list(map(lambda a,d: SquareVectorNorm(np.cross(sol - a, d)), a, d))
        dists = VectorNorm(dists)*np.sqrt(1/len(p))
        dists = dists.item()  # Convert to regular float
        return([sol.tolist(),dists])

# 3D dimensional voxel traversal, gives cell-indices back starting from point p and moving in v direction, subject two cell-bounds bounds
def DirectionalVoxelTraversal2(p,v,bounds,logfile=''):
    if len(p) == len(v) == len(bounds) == 3:    # p, v, and bounds should be 3 dimensional
        curindex = list(map(PositionSorted,bounds,p))
        if (-1 in curindex or VectorNorm(v)==0):
            if(logfile != ''):
                flog = open(logfile, 'a')
                flog.write("ray starts outside bounds!\np:" + str(p) + "\nv:" + str(v) + "\ncell index:" + str(curindex) + "\n")
                flog.close()
            print("ray starts outside bounds!")
            print("p:",p)
            print("v:",v)
            print("cell index:",curindex)
            return([])
        else:
            direction = list(map(Sgn,list(v)))
            relbounds = list(map(lambda a,b:list([i-b for i in a]),list(bounds),list(p)))
            times = list(map(lambda a,b:list([SpecialDivision(i,b) for i in a]),relbounds,v))
            times = list(map(lambda a,b:list([[i,b,False] for i in a]),times,[0,1,2]))
            for i in range(3):
                if(direction[i]==1):
                    times[i][-1][2] = True   
                else:
                    times[i][0][2] = True   
            times = joinlists(times)
            times = list(filter(lambda x: x[0]>0,times))  # Should be > !
            times = list(sorted(times,key=lambda x: x[0]))
            times = list(it.takewhile(lambda x: x[2]==False,times))
            times = list([x[1] for x in times])
            out = [copy.copy(curindex)]
            for index in times:
                new = list(out[-1]);
                new[index] += direction[index]
                out.append(list(new)) 
            return(out)
    else:
        print("dimension mismatch!")
        if(logfile != ''):
            flog = open(logfile, 'a')
            flog.write("dimension mismatch!\n")
            flog.close()
        return([])

# 3D dimensional voxel traversal, gives cell-indices back starting from point p and moving in v direction, subject two cell-bounds bounds
def DirectionalVoxelTraversal(p,v,bounds,logfile=''):
    if len(p) == len(v) == len(bounds) == 3:    # p, v, and bounds should be 3 dimensional
        curpos = list(p)
        cellindex = list(map(PositionSorted,bounds,p))
        if (-1 in cellindex or VectorNorm(v)==0):
            if(logfile != ''):
                flog = open(logfile, 'a')
                flog.write("ray starts outside bounds!\np:" + str(p) + "\nv:" + str(v) + "\ncell index:" + str(cellindex) + "\n")
                flog.close()
            print("ray starts outside bounds!")
            print("p:",p)
            print("v:",v)
            print("cell index:",cellindex)
            return([])
        else:
            sgns = list(map(Sgn,v))
            sgnspart = list(map(lambda x: 1 if x == -1 else 0, sgns))
            
            cont = True
            steps = 0;
            out = [copy.copy(cellindex)]
            while(cont and steps < 10000): # Steps is just a safety thing for now, can be removed later
                steps += 1
                newindices = list(map(lambda x,y: x+y, cellindex, sgns))
                newbounds = [0,0,0]
                for i in range(3):
                    if (0 <= newindices[i]  <= len(bounds[i])-1):  # Should there be + sgnspart[i] in the middle term?
                        newbounds[i] = bounds[i][newindices[i] + sgnspart[i]]
                    else:
                        newbounds[i] = math.inf  # Bounds are at infinity, so the 'next' bounds are infinitely far away
            
                if (math.inf not in newbounds):
                    ts = [0,0,0]
                    for i in range(3):
                        if (v[i]==0):
                            ts[i] = math.inf     # It will take infinite amount of time
                        else:
                            ts[i] = (newbounds[i]-curpos[i])/v[i]
                
                    order = sorted(range(len(ts)), key=lambda k: ts[k])  # Find the 'times' needed to the next boundaries in each dimensions
                    minpos = order[0]                                    # Find the dimension for which the time is shortest
                    
                    cellindex[minpos] += sgns[minpos]
                    ts = ts[minpos]
                    for i in range(3):
                        curpos[i] += v[i]*ts
        
                    out.append(copy.copy(cellindex))
                else:
                    #print("at edge")
                    cont = False
        return(out)
    else:
        print("dimension mismatch!")
        if(logfile != ''):
            flog = open(logfile, 'a')
            flog.write("dimension mismatch!\n")
            flog.close()
        return([])

def AtFace(bmin, bmax, hitb):
    return(bmin <= hitb <= bmax)

# Projects ray (defined by p,v) on to an AABB (axis aligned bounding box).
# [boolhit, boolinside, pnew, v] boolhit tells if it hits AABB, abd boolinside tells if it is projected, and pnew the new position or [] in case it misses.
# Note that ray can be projected on to an AABB with negative 'time'...
def PrepareRay(p,v,bounds):
    xmin = bounds[0][0]
    xmax = bounds[0][1]
    ymin = bounds[1][0]
    ymax = bounds[1][1]
    zmin = bounds[2][0]
    zmax = bounds[2][1]
    x = p[0]
    y = p[1]
    z = p[2]
    newv = normalize(v).tolist()       # v is normalized
    vx, vy, vz = newv
    if(xmin < x < xmax and ymin < y < ymax and zmin < z < zmax):
        return([True, True, p, newv])  # Return False and original point
    else:
        t = list(map(SpecialDivision,[xmin - x, xmax - x, ymin - y, ymax - y, zmin - z, zmax - z],[vx, vx, vy, vy, vz, vz]))
        ip = [0 for tmp in range(6)]
        for i in range(6):
            ti = t[i]
            if(abs(ti) == math.inf):
                ip[i] = [math.inf, math.inf, math.inf]
            else:
                ip[i] = [x + vx*ti, y + vy*ti, z + vz*ti]

    atfacex1 = AtFace(ymin, ymax, ip[0][1]) and AtFace(zmin, zmax, ip[0][2])
    atfacex2 = AtFace(ymin, ymax, ip[1][1]) and AtFace(zmin, zmax, ip[1][2])
    atfacey1 = AtFace(xmin, xmax, ip[2][0]) and AtFace(zmin, zmax, ip[2][2])
    atfacey2 = AtFace(xmin, xmax, ip[3][0]) and AtFace(zmin, zmax, ip[3][2])
    atfacez1 = AtFace(xmin, xmax, ip[4][0]) and AtFace(ymin, ymax, ip[4][1])
    atfacez2 = AtFace(xmin, xmax, ip[5][0]) and AtFace(ymin, ymax, ip[5][1])
    data = [t, [atfacex1, atfacex2, atfacey1, atfacey2, atfacez1, atfacez2], ip]
    data = list(zip(*data))                                   # Data will be a list of lists, each having the form: time-till-hit, hits face (boolean), position it hits a plane
    data = list(filter(lambda xx: xx[1]==True, data))         # Select only those that hit a face #don't change to 'is True' numpy bool possibility
    if(len(data)>0):
        data = sorted(data, key = lambda x: x[0])             # Sort by arrival time (time till hit)
        return([True, False, data[0][2], newv])               # Position it hits the plane of first-hit
    else:
        return([False, False, [], newv])
    
    
    
def SpaceTraversalMatching(raydata, boundingbox, nx = 75, ny = 75, nz = 75, cammatchfunc = lambda x: len(x)>2, maxmatchesperray = 2, maxdistance = 999.9, neighbours = 6, logfile = ''):

    if(len(boundingbox)==3 and len(boundingbox[0]) == 2 and len(boundingbox[1]) == 2 and len(boundingbox[2]) == 2 and nx >= 5 and ny >= 5 and nz >= 5):
        
        if(neighbours == 0):
            neifhbours = [[0,0,0]]
        elif(neighbours == 6):
            neighbours = [[-1,0,0],[0,-1,0],[0,0,-1],[0,0,1],[0,1,0],[1,0,0],[0,0,0]]
        elif(neighbours == 18):
            neighbours = [[-1,-1,0],[-1,0,-1],[-1,0,0],[-1,0,1],[-1,1,0],[0,-1,-1],[0,-1,0],[0,-1,1],[0,0,-1],[0,0,1],[0,1,-1],[0,1,0],[0,1,1],[1,-1,0],[1,0,-1],[1,0,0],[1,0,1],[1,1,0],[0,0,0]]
        elif(neighbours == 26):
            neighbours = [[-1,-1,-1],[-1,-1,0],[-1,-1,1],[-1,0,-1],[-1,0,0],[-1,0,1],[-1,1,-1],[-1,1,0],[-1,1,1],[0,-1,-1],[0,-1,0],[0,-1,1],[0,0,-1],[0,0,1],[0,1,-1],[0,1,0],[0,1,1],[1,-1,-1],[1,-1,0],[1,-1,1],[1,0,-1],[1,0,0],[1,0,1],[1,1,-1],[1,1,0],[1,1,1],[0,0,0]]
        if(type(neighbours) == list):
            if(logfile != ''):
                def LOGprint(*out):
                    flog = open(logfile, 'a')
                    flog.write(" ".join(map(str,list(out))) + "\n")
                    flog.close()
            else:
                def LOGprint(*out):
                    print(*out)

            LOGprint(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
            bounds = list(map(lambda x, n: list([x[0]+i*(x[1]-x[0])/n for i in range(n+1)]),boundingbox,[nx,ny,nz]))
            for i in range(3):
                bounds[i][-1] = boundingbox[i][-1]
            
            LOGprint("# of cells:",nx*ny*nz)
            rays = list(raydata)

            # Prepare the rays so that they are inside the box or on the side.
            cammarkerfunc = lambda x: x[0]        # First element is camera ID
            raymarkerfunc = lambda x: x[1]        # Second element is the ray ID
            raydb = {}                            # Store a dictionary of (cameraid, rayid): [pos, direction]
            validrays = []                        # Store the transformed rays                 
            numrays = col.Counter()               # Store the number of rays per camera
            invalidcounter = col.Counter()        # Store the number that are invalid
            for r in rays:
                camid = cammarkerfunc(r)
                rayid = raymarkerfunc(r)
                numrays[camid] += 1

                pp = list(r[2:5])
                vv = list(r[5:8])
                out = PrepareRay(pp, vv, boundingbox)

                if(out[0]):  # If it does not miss (hit or inside)
                    raydb[(camid, rayid)] = [out[2], out[3]]                     # Tuple of cam-ray ids
                    validrays.append([camid, rayid, out[1], out[2], out[3]])     # CamID, rayID, bool inside, pnew, v 
                else:
                    #raydb[(camid, rayid)] = [copy.copy(pp), copy.copy(vv),"not used"]      # Store in raydb even if it misses the bounding box, mark them
                    invalidcounter[camid] += 1

            #for k, v in raydb.items():
            #    print(VectorNorm(v[1]))
            
            
            LOGprint("# of rays for each camera:", dict(numrays),"\n# of rays that miss the bounding box:", dict(invalidcounter))
            if(len(dict(numrays))>10):
                LOGprint("# of cameras is large:", len(dict(numrays)) ," Double check the input!")
                return([])
                         

            traversed = []
            for ray in validrays:
                if(ray[2]):  # Ray is inside, traverse both forward and backward
                    out = DirectionalVoxelTraversal2(list(ray[3]),list(ray[4]),bounds,logfile) + \
                          DirectionalVoxelTraversal2(list(ray[3]),list(map(lambda x: -1*x, ray[4])),bounds,logfile)
                else:        # Ray is at edge, traverse in forward direction only
                    out = DirectionalVoxelTraversal2(list(ray[3]),list(ray[4]),bounds,logfile)

                out = uniqify(expandallneighbours(out,neighbours),lambda x: tuple(x))    # Expand in neighbourhood and remove duplicates
                ext = [[ray[0],ray[1],o] for o in out]                                   # Creates long list of camid, rayid, cellindex
                traversed.extend(ext)                                                    # Pile up these into traversed

            LOGprint("# of voxels traversed after expansion:", len(traversed))
            cellfunc = lambda x: x[2]
            traversed = list(sorted(traversed,key=cellfunc))                             # Sort based on cell. Sort is needed before groupby
            traversed = [list(g) for k, g in it.groupby(traversed, cellfunc)]            # Group elements by same cell
            LOGprint("Sorted and grouped by cell index. # of groups:",len(traversed))
            traversed = list(filter(cammatchfunc,traversed))                             # Prune based on number of rays (fast rough filter, cam filter later)
            LOGprint("Rough pruned based on number of cameras:",len(traversed))
            traversed = maplevel(lambda x:[x[0],x[1]], traversed, 2)                     # Remove cellindex, not needed anymore, leave [camid, rayid]
            traversed = list(map(lambda x:[list(g) for k, g in it.groupby(x, cammarkerfunc)],traversed))
            #LOGprint("Cell index removed and grouped by camera for each cell")
            traversed = list(filter(cammatchfunc,traversed))                             # Prune based on number of different cameras
            LOGprint("Pruned based on number of cameras:",len(traversed))
            candidates = copy.copy(list(map(lambda x:list(list(tup) for tup in it.product(*x)),traversed)))        # All combinations between all cameras
            candidates = joinlists(candidates)                                           # Flatten a list of lists to a single list
            LOGprint("Flattened list of candidates:",len(candidates))
            candidates = uniqify(candidates,lambda x: tuple(list(joinlists(x))))         # Delete duplicates, flattened list as tag
            LOGprint("Duplicate candidates removed:",len(candidates))
            candidates = sorted(candidates)

                        
            #LOGprint("Computing match position and quality of candidates...")
            newcandidates = []
            for c in candidates:
                pvdata = list([raydb[tuple(x)] for x in c])
                pdata = list([x[0] for x in pvdata])
                vdata = list([x[1] for x in pvdata])   
                out = ClosestPointToLines(pdata,vdata)
                newcandidates.append([c,list(out[0]),out[1]])
                
            #hullpts = [[20.0,5.0,165.0],[20.0,10.0,160.0],[20.0,10.0,165.0],[20.0,15.0,160.0],[20.0,15.0,165.0],[20.0,20.0,160.0],[20.0,20.0,165.0],[20.0,25.0,160.0],[20.0,25.0,165.0],[25.0,0.0,165.0],[25.0,0.0,170.0],[25.0,5.0,155.0],[25.0,10.0,155.0],[25.0,15.0,155.0],[25.0,20.0,155.0],[25.0,25.0,155.0],[25.0,25.0,175.0],[25.0,30.0,170.0],[25.0,35.0,155.0],[25.0,35.0,160.0],[30.0,-5.0,170.0],[30.0,0.0,160.0],[30.0,5.0,150.0],[30.0,10.0,150.0],[30.0,15.0,150.0],[30.0,20.0,150.0],[30.0,25.0,150.0],[30.0,25.0,180.0],[30.0,30.0,150.0],[30.0,30.0,180.0],[30.0,35.0,175.0],[30.0,40.0,165.0],[30.0,45.0,155.0],[35.0,-10.0,175.0],[35.0,-5.0,165.0],[35.0,-5.0,180.0],[35.0,0.0,155.0],[35.0,10.0,145.0],[35.0,15.0,145.0],[35.0,20.0,145.0],[35.0,25.0,185.0],[35.0,30.0,185.0],[35.0,35.0,150.0],[35.0,35.0,185.0],[35.0,40.0,180.0],[35.0,45.0,155.0],[35.0,45.0,170.0],[35.0,50.0,160.0],[40.0,-10.0,175.0],[40.0,-10.0,180.0],[40.0,-5.0,165.0],[40.0,0.0,155.0],[40.0,5.0,145.0],[40.0,10.0,145.0],[40.0,15.0,145.0],[40.0,30.0,150.0],[40.0,40.0,180.0],[40.0,45.0,155.0],[40.0,45.0,175.0],[40.0,50.0,160.0],[40.0,50.0,165.0],[45.0,-5.0,175.0],[45.0,0.0,165.0],[45.0,5.0,155.0],[45.0,10.0,150.0],[45.0,15.0,150.0],[45.0,40.0,180.0],[45.0,50.0,160.0],[45.0,50.0,165.0],[50.0,15.0,160.0],[50.0,20.0,160.0],[50.0,25.0,175.0],[50.0,30.0,175.0],[50.0,35.0,175.0],[50.0,45.0,165.0],[55.0,15.0,170.0],[55.0,20.0,170.0],[55.0,25.0,170.0],[55.0,30.0,170.0],[55.0,35.0,170.0],[55.0,40.0,170.0]];
            #delaun = sps.Delaunay(hullpts)  # Define the Delaunay triangulation            
            #inq = delaun.find_simplex([x[1] for x in newcandidates])>0
            #inq = inq.tolist();
            #print(inq)
            #print(type(inq))
            #newcandidates = list(map(lambda x,y: x + [y],newcandidates,inq))

            candidates = copy.copy(newcandidates);
            del newcandidates
            #LOGprint("Sorting candidate matches by quality of match...")
            candidates = sorted(candidates,key=lambda x: (-len(x[0]), x[2]))       #sort by number of cameras then by error
            #LOGprint("Here are upto 9999 of the best matches:")
            #LOGprint("Index, [camid rayid ....] Position, Mean square distance")
            #print("num candidates:",len(candidates))
            #for i in range(1,len(candidates),100):
            #    print(i,candidates[i])

            LOGprint("Selecting the best matches with upto",maxmatchesperray,"match(es)/ray out of",len(candidates),"candidates")
            #now we want to pick the best matches first and match each ray at most maxmatchesperray
            approvedmatches = []             # Store approved candidates
            matchcounter = col.Counter()     # Keep track of how many they are matched
            for cand in candidates:
                if(cand[2] < maxdistance):
                    valid = True;
                    for idpair in cand[0]:
                        if(matchcounter[tuple(idpair)]>=maxmatchesperray):
                            valid = False
                            #print(idpair,"has been matched already",maxmatchesperray,"time(s)")
                            break
                    if(valid==True):
                        for idpair in cand[0]:
                            matchcounter[tuple(idpair)] += 1

                        approvedmatches.append(list(cand))

            LOGprint("Selecting done.",len(approvedmatches),"matched found (out of",len(candidates),"candidates)")
            #print("Here are the approved matches:")
            #print("Index, [camid rayid ....] Position, Mean square distance")    
            
            return(approvedmatches)
        else:
            print("Neighbours should be a 0, 6, 18, or 26 or a list of triplets:",neighbours)  
            return([])
    else:
        print("Something went wrong:")
        print("Bounding box should be of the form [[xmin,xmax],[ymin,ymax],[zmin,zmax]]",boundingbox)
        print("nx,ny,nz should be >=5 otherwise you will do a lot of comparisons!",[nx,ny,nz])    
        return([])
    
