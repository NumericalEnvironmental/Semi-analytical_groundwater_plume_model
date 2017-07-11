##################################################################################
#
# 2-D analytical advective-dispersive solute transport
#
# - basis is instantaneous point source model
# - integrated in time to represent time-dependent mass influx
# - additional integration over finite source area spatial footprint
#
# by Walt McNab
#
##################################################################################



### data types ###


type Aquifer
	b::Float64 						# aquifer thickness
	phi::Float64 					# porosity
	v::Float64 						# mean pore velocity
	theta::Float64 					# groundwater flow direction, counterclockwise from due east (radians)
	R::Float64 						# retardation coefficient
	k::Float64 						# solute first-order degradation rate constant		
	fl::Float64 					# ratio of longitudinal dispersivity to plume length scale (e.g., 0.1)
	ft::Float64 					# ratio of transverse dispersivity to plume length scale
	Dlmin::Float64 					# minimum longitudinal dispersivity
	Dtmin::Float64 					# minimum transverse dispersivity	
end


type Source
	name::AbstractString 			# name of source area, used to designate time series input file
    x0::Float64 					# center point location (absolute coordinate); corresponds to (0, 0) in relative coordinates
    y0::Float64
    tend::Float64					# time at (abrupt) cessation of source flux
    xs::Float64 					# start of source area w.r.t. flow direction (relative coordinate)
    xf::Float64						# end of source area w.r.t. flow direction (relative coordinate)
end


### transport modeling functions ###


function Cslug(x::Float64, y::Float64, t::Float64, v::Float64, M::Float64, phi::Float64, b::Float64, Dl::Float64, Dt::Float64, k::Float64)::Float64
	# Bear, 1972 solute transport analytical solution
    C = M/(4.*pi*phi*b*sqrt(Dl*Dt)*t) * exp(-((x-v*t)^2)/(4.*Dl*t) - y^2/(4.*Dt*t) - k*t)
    return C
end


function dCdt(tau::Float64, x::Float64, y::Float64, teval::Float64, aquifer::Aquifer)::Float64
	# slice-in-time for a small [point] source
	# tau = integration value (along source history; absolute value for time)
	# teval = elapsed time since time zero
	# x and y are relative offsets from source location
	t = teval - tau 				# lag time in evaluation time versus source term
	Dl, Dt = Dispers(tau, aquifer.v, aquifer.R, aquifer.fl, aquifer.ft, aquifer.Dlmin, aquifer.Dtmin)
	dM = tfunc(t)
	return Cslug(x, y, tau, aquifer.v/aquifer.R, dM/aquifer.R, aquifer.phi, aquifer.b, Dl, Dt, aquifer.k)
end
	
	
function Ct(xs::Float64, ys::Float64, xp::Float64, yp::Float64, t::Float64, aquifer::Aquifer, source::Source)::Float64
	# concentration from a point source stemming from time-integration
	# xs and ys are source locations (relative coordinates)
	# xp and yp are monitoring point locations (relative coordinates)
	return quadgk(tau -> dCdt(tau, xp-xs, yp-ys, t, aquifer), t-source.tend, t)[1]
end

	
function Cy(x::Float64, xp::Float64, yp::Float64, t::Float64, aquifer::Aquifer, source::Source)::Float64	
	# x is the position along the source zone x-line (relative coordinates), supplied as an integration variable
	# y is the variable of integration for a line running between xbotf and xtopf (relative coordinates), as defined by f's(x)
	# xp and yp are monitoring point locations (relative coordinates)
	ys = xbotf(x)
	yf = xtopf(x)
	return quadgk(y -> Ct(x, y, xp, yp, t, aquifer, source), ys, yf)[1]
end
	
	
function C(xp::Float64, yp::Float64, t::Float64, aquifer::Aquifer, source::Source)::Float64
	# compute integrated concentration at location (x, y) and time t for this particular source
	return quadgk(x -> Cy(x, xp, yp, t, aquifer, source), source.xs, source.xf)[1]
end   


function Dispers(t::Float64, v::Float64, R::Float64, fl::Float64, ft::Float64, Dlmin::Float64, Dtmin::Float64)
    # compute dispersion coefficient component as a function of plume length scale
	Dl = fl * v^2 * t / R + Dlmin
	Dt = ft * v^2 * t / R + Dtmin 		
	return Dl, Dt
end



### utility functions (input/output, etc.) ###


function GetAquifer()::Aquifer
	# read/interpret aquifer properties from file
	data = readdlm("aquifer.txt", '\t', header=false)
    b = Float64(data[1, 2]) 				# aquifer thickness
    phi = Float64(data[2, 2]) 				# porosity
    v = Float64(data[3, 2]) 				# mean pore velocity
    theta = Float64(data[4, 2]) 			# groundwater flow direction
    R = Float64(data[5, 2]) 				# retardation coefficient
    k = Float64(data[6, 2]) 				# solute first-order degradation rate constant	
    fl = Float64(data[7, 2]) 				# ratio of longitudinal dispersivity to plume length scale (e.g., 0.1)
    ft = Float64(data[8, 2]) 				# ratio of transverse dispersivity to plume length scale
	aLmin = Float64(data[9, 2]) 			# minimum longitudinal dispersivity
	aTmin = Float64(data[10, 2]) 			# minimum transverse dispersivity	
	aquifer = Aquifer(b, phi, v, theta, R, k, fl, ft, aLmin*v, aTmin*v)
	println("Read aquifer properties.")
	return aquifer
end


function GetSources()::Array{Source, 1}
	# read and process source term parameters from file
	sources = []
	data = readdlm("sources.txt", '\t', header=true)
	for i = 1:size(data[1], 1) 							# can have multiple sources
		name = data[1][i, 1]
        x0 = Float64(data[1][i, 2])
        y0 = Float64(data[1][i, 3])
        tend = Float64(data[1][i, 4])
        xs = Float64(data[1][i, 5])
        xf = Float64(data[1][i, 6])
		push!(sources, Source(name, x0, y0, tend, xs, xf))
	end
	println("Read source term properties.")
	return sources	
end


function GetGrid()
	# read and process grid constraints from file; return set of points
	xgrid = []
	ygrid = []
	data = readdlm("grid.txt", '\t', header=true)
    x0 = Float64(data[1][1, 2])
    y0 = Float64(data[1][1, 3])
    xf = Float64(data[1][2, 2])
    yf = Float64(data[1][2, 3])
    numx = Int64(data[1][3, 2])        
    numy = Int64(data[1][3, 3])  
	dx = (xf-x0)/numx
	dy = (yf-y0)/numy
	for i = 1:numx + 1, j = 1:numy + 1
		push!(xgrid, x0 + (i-1)*dx)
		push!(ygrid, y0 + (j-1)*dy)		
	end	
	return xgrid, ygrid 
end


function WriteOutput(xgrid, ygrid, sources, concs)
	# write out modeled concentration values for grid points
	fname = "concs_output.csv"
	csvfile = open(fname,"w")
	line_out = "x" * "," * "y"
	for src in sources
		line_out *= "," * src.name
	end
	line_out *= "," * "total"
	println(csvfile, line_out)
	for i = 1:length(xgrid)
		line_out = string(xgrid[i]) * "," * string(ygrid[i])
		for j = 1:length(sources)
			line_out *= "," * string(concs[i,j])
		end
		line_out *= "," * string(sum(concs[i,:]))
		println(csvfile, line_out)		
	end
	close(csvfile)		
	println("Wrote output.")
end


function Rotate(x, y, theta)
    xprime = x*cos(theta) - y*sin(theta)
    yprime = x*sin(theta) + y*cos(theta)		
    return xprime, yprime
end


### read plume source terms; this must be implemented as global scope ###


tfunc_string = []
xtopf_string = []
xbotf_string = []
data = readdlm("sources.txt", '\t', header=true)
for i = 1:size(data[1], 1) 							# can have multiple sources
	push!(tfunc_string, data[1][i, 7])
	push!(xtopf_string, data[1][i, 8])
	push!(xbotf_string, data[1][i, 9])		
end


### main script ###


function ShapeSource(time_output)

	# read in model parameter sets
	aquifer = GetAquifer() 						# read aquifer properties
	sources = GetSources() 						# read source term properties
	xgrid, ygrid = GetGrid() 					# read grid settings
	
	# calculate
	concs = zeros(Float64, length(xgrid), length(sources))
	for (j, src) in enumerate(sources)
	println("Source = " * "\t" * src.name)
		# redefine global functions to match sources
		@eval tfunc(t) = $(parse(tfunc_string[j]))  	
		@eval xtopf(x) = $(parse(xtopf_string[j]))
		@eval xbotf(x) = $(parse(xbotf_string[j]))
		for i = 1:length(xgrid)
			xprime, yprime = Rotate(xgrid[i]-src.x0, ygrid[i]-src.y0, -aquifer.theta)
			concs[i, j] = C(xprime, yprime, time_output, aquifer, src)
		end
	end
	
	# write to output file
	WriteOutput(xgrid, ygrid, sources, concs)
	
	println("Done.")
	
end


### run script ###

time_output = 50.1
@time ShapeSource(time_output)
