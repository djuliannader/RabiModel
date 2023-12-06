push!(LOAD_PATH, pwd())
module reading
export stringtofloatlist

function stringtofloatlist(x::String)
	 comas=Vector{Int64}()
	 for i in 1:length(x)
	     if x[[i]]==","
	     push!(comas,i)
	     end
	 end
	 if length(comas)==0
	    lista=[parse(Float64,x[2:length(x)-1])]
	    return lista
	 end
         lista=Vector{Float64}()
	 ni=2
	 nf=comas[1]
	 nl=parse(Float64,x[ni:nf-1])
	 append!(lista,nl)
	 if length(comas)>1
	    for i in 1:length(comas)-1
	       ni=comas[i]
	       nf=comas[i+1]
	       nl=parse(Float64,x[ni+1:nf-1])
	       append!(lista,nl)
	       end
         end
	 ni=nf
	 nf=length(x)-1
	 nl=parse(Float64,x[ni+1:nf])
	 append!(lista,nl)
	 return lista
	 end

#y="[3.666678,2.563,3.36699,6.38,8.45]"
#s=readinglist(y)
#println(s[3])

end