using Comodo
using Comodo.DelaunayTriangulation
using GLMakie
using GeometryBasics
using FileIO
using Statistics
using LinearAlgebra
using Geogram

# 
Z_cut_offset_top = 200.0
Z_cut_level_skin = 160.0 #(100-160)
Z_cut_level_tibia = Z_cut_level_skin+50
Z_cut_level_fibula = Z_cut_level_skin+10

Z_thickness_distal = 20 #distal thickness can be changed

pointSpacing = 5.0 
fileName_set = (
     "/home/simulimb/Che_simamp/SimulateAmputation.jl/Che_assets/FMA7163_right_leg_isolated_remesh_25k.stl",
"/home/simulimb/Che_simamp/SimulateAmputation.jl/Che_assets/FMA24474.stl",
"/home/simulimb/Che_simamp/SimulateAmputation.jl/Che_assets/FMA24486.stl",
"/home/simulimb/Che_simamp/SimulateAmputation.jl/Che_assets/FMA24477.stl",
"/home/simulimb/Che_simamp/SimulateAmputation.jl/Che_assets/FMA24480.stl")

nameSet = ("skin","femur","patella","tibia","fibula")
indSkin = 1
indFemur = 2
indPatella = 3
indTibia = 4
indFibula = 5

# nameSetCut_top = ("skin")

nameSetCut_top = ("skin","femur")
nameSetCut_bottom = ("skin","tibia","fibula")
nameSetCut_smooth = ("femur","tibia","fibula")

n = normalizevector(Vec{3, Float64}(0.0,0.0,1.0))# Cutting plane normal
snapTolerance = 1e-6

# FMAID 	English name
# 7163  	Skin
# 16586 	Right hip bone
# 24474 	Right femur
# 24486     Right patella
# 24477 	Right tibia
# 24480 	Right fibula

FF = []
VV = []
for i in 1:1:length(fileName_set)
    M = load(fileName_set[i])
    F = tofaces(faces(M))
    V = topoints(coordinates(M))
    V = [Point{3,Float64}(v) for v in V]
    
    F,V,_,_ = mergevertices(F,V; pointSpacing=pointSpacing)

    # Remesh evenly using desired point spacing

    # Determine number of points to request to closely match desired point spacing
    np = spacing2numvertices(F,V,pointSpacing) 

    # Remesh with desired number of points 
    F,V = ggremesh(F,V; nb_pts=np,suppress = true, cleanup = false)

    push!(FF,F)
    push!(VV,V)
end

FF_ori = deepcopy(FF)
VV_ori = deepcopy(VV)

# Cutting surfaces 
#finding the midpoint of patella 
v_mid_patella = mean(VV[indPatella])
cut_plane_top = v_mid_patella[3]+Z_cut_offset_top
cut_plane_bottom = v_mid_patella[3]-Z_cut_level_skin
cut_plane_tibia = v_mid_patella[3] -Z_cut_level_tibia
cut_plane_fibula = v_mid_patella[3] -Z_cut_level_fibula
p = [v_mid_patella[1],v_mid_patella[2],v_mid_patella[3]+Z_cut_offset_top] # Point on cutting plane for the top skin 
for modelName in nameSetCut_top
    indexNow = findall(nameSet.==modelName)[1]
    Fn = FF[indexNow]
    Vn = VV[indexNow]
   B = [all([Vn[i][3] < cut_plane_top for i in f]) for f in Fn]
    B_new = fixteeth(Fn,B,method=:add)       
    Fn1 = Fn[B_new]
    Fn1,Vn1,_ = remove_unused_vertices(Fn1,Vn)
    FF[indexNow] = Fn1
    VV[indexNow] = Vn1
end
#skin bottom cut 
cut_level = [cut_plane_bottom, cut_plane_tibia,cut_plane_fibula]
for modelName in nameSetCut_bottom
    indexNow = findall(nameSet.==modelName)[1]
    i = findall(nameSetCut_bottom .== modelName)[1]
    Fn = FF[indexNow] 
    Vn = VV[indexNow]  
    B = [all([Vn[j][3] > cut_level[i] for j in f]) for f in Fn]
    B_new = fixteeth(Fn,B,method=:remove)       
    Fn1 = Fn[B_new]
    Fn1,Vn1,_ = remove_unused_vertices(Fn1,Vn)
    FF[indexNow] = Fn1
    VV[indexNow] = Vn1
end 
#end tibia 
V_tibia = VV[indTibia]
Eb_tibia = boundaryedges(FF[indTibia])
indCurve_tibia = edges2curve(Eb_tibia)
indCurve_tibia = indCurve_tibia[1:end-1]
v_mid_low_tibia = mean(V_tibia[indCurve_tibia])
p_distal = Point{3,Float64}(v_mid_low_tibia[1],v_mid_low_tibia[2],v_mid_low_tibia[3]-Z_thickness_distal)

#zlevel top 
zlevel = v_mid_patella[3]

#smoothen the skin
#isolate the skin 
Fs = FF[1]
Vs = VV[1]
# #find the boundary edges (top and bottom )
Eb = boundaryedges(Fs)
# indCurve_2 = edges2curve(Eb)
# #extract the bottom edges: 
# #the third coordinate is compared and lower than a certain threshold is removed  
Eb_low = Vector{eltype(Eb)}()
Eb_high = Vector{eltype(Eb)}()
for e in Eb         
    if Vs[e[1]][3] < zlevel        
      push!(Eb_low,e)
    else push!(Eb_high,e)  
    end
end
indCurve = edges2curve(Eb_low)
indCurve = indCurve[1:end-1] 

#top surface:

indCurve_1 = edges2curve(Eb_high)
Vtop = Vs[indCurve_1[1:end-1]]

meanZ = mean([v[3] for v in Vtop])
Z = [v[3] for v in Vtop]
maxZ = maximum(Z)
maxid_z =findall(z->z==maxZ,Z)
minZ = minimum(Z)
minid_z =findall(z->z==minZ,Z)
limit = maxZ- ((maxZ - minZ)/pointSpacing)
# indRemove = findall(z->z > meanZ+1e-3,Z)
B = [all([Vs[i][3] <= meanZ+.1 for i in f]) for f in Fs]
B_new = fixteeth(Fs,B,method=:add)       
Fs1 = Fs[B_new]
Fs1,Vs1,_ = remove_unused_vertices(Fs1,Vs)
    


ids = findall(z->z> limit && z < maxZ,Z)
Vtop_shift = [Point{3,Float64}(v[1], v[2], meanZ) for v in Vtop]
#Visualization
Vs[indCurve_1[1:end-1]] .= Vtop_shift
# cons_seg = deepcopy(indCurve_1)
# cons_seg = reverse!(cons_seg)
# TRn = triangulate(Vtop; boundary_nodes=cons_seg, delete_ghosts=true)
# # Fs1 = [TriangleFace{Int}(tr) for tr in each_solid_triangle(TRn)] 
# Vtop1 = get_points(TRn)
        

fig1 = Figure(size=(1200,800))
ax1 = LScene(fig1[1,1])
for i in 1:1:length(fileName_set)
    hp1 = poly!(ax1,GeometryBasics.Mesh(VV[i],FF[i]), color=:white, shading = FastShading, transparency=false,strokecolor=:black,strokewidth=1)
    
end
hp2 = wireframe!(ax1,GeometryBasics.Mesh(Vs,Eb_high), linewidth=5,color=:red)
# hp3 = scatter!(ax1,Vtop[maxid_z], markersize=20,color=:blue)
# hp4 = scatter!(ax1,Vtop[minid_z], markersize=20,color=:green)
# hp5 = scatter!(ax1,Vtop[indRemove], markersize=20,color=:yellow)
ax2 = LScene(fig1[1,2])
hp1 = poly!(ax2,GeometryBasics.Mesh(Vs1,Fs1), color=:white, shading = FastShading, transparency=false,strokecolor=:black,strokewidth=1)
hp2 = wireframe!(ax2,GeometryBasics.Mesh(Vs,Eb_high), linewidth=5,color=:red)

hp6 = scatter!(ax2,Vtop, markersize=20,color=:black)
fig1