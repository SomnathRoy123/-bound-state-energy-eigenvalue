using PyPlot
function secant(x1,x2,acc,f,g,k)
    iter=0 #set iteration numebr to zero
    #check if roots are at initial guesses
    if f(g,x1,k)==0
        x1
    elseif f(g,x2,k)==0
        x2
    #if not then start the iteration process of secant method
    else
        guess=x1 #set initial guess to x1
        while(abs(f(g,guess,k))>acc && iter<100)
            guess=x2-f(g,x2,k)*((x2-x1)/(f(g,x2,k)-f(g,x1,k))) #keep updating guess using the secant rule formula
            x1=x2 #shift guess x2 to x1 and new guess to x2 and use them again to find new geuss
            x2=guess
            iter=iter+1
            #update iteration number
        end
        guess #return the guess
    end
end
function bisect(lo,hi,acc,f,e) #usual bisection function
    iter=0
    if f(lo,e)==0
        lo
    elseif f(hi,e)==0
        hi
    elseif f(lo,e)*f(hi,e)>0
        print("Interval is not suitable for bisection method!!")
    else
        guess=(hi+lo)/2
        while (iter<100)
            guess=(lo+hi)/2
            iter=iter+1
            if f(lo,e)*f(guess,e)<0
                hi=guess
            else
                lo=guess
            end
            if (abs(hi-lo)<acc)
                break
            end
        end
        guess
    end
end
root=2^(1/6)
acc=1.0e-16
#lennard jonnes function to find root of v(x)=e 
function LJ(x,e)
    4(1/x^12 - 1/x^6)-e
end
#function to find the integral in order to find s(e)
function LJp(x,e)
    sqrt(abs(e-4(1/x^12 - 1/x^6)))
end
#functions that find r_in and r_out
function r_infind(e)
    (bisect(1,root,acc,LJ,e))
end
function r_outfind(e)
    (bisect(root,10,acc,LJ,e))
end
#simpson numerical integrator function
function simp(f,lo,hi,e)
    xs=range(lo,stop=hi,step=(hi-lo)/20000)
    ys=f.(xs,e)
    (4*sum(ys[2:2:end-1])+2*sum(ys[3:2:end-2])+(ys[1]+ys[end]))*((hi-lo)/60000)
end
#s(e) defined in an analytic form with variable gamma (g) and k which takes values π,3π,5π...
function s(g,e,k)
    g*simp(LJp,r_infind(e),r_outfind(e),e)-k
end
i=0
#vector to store eigenvalues for H2
H2=Vector{Float64}()
print("Hydrogen:\n")
print("Ground state: \n")
#manually finding the two states for H2
m=secant(-0.96,-0.94,1e-10,s,21.7,(i+0.5)*2*π)
push!(H2,m)
print(m,"\t value of s(e)= ",s(21.7,m,0),"\n\n")
print("First excited state: \n")
m=secant(-0.15,-0.12,1e-10,s,21.7,(1+0.5)*2*π)
push!(H2,m)
print(m,"\t value of s(e)= ",s(21.7,m,0),"\n\n")
print("No more states for hydrogen could be found using this program!\n")
print("Oxygen:\n")
#vector to store the states of O2
O2=Vector{Float64}()
#loop to compute states of O2
while(i<6)
    a=-0.96
    b=-0.94
    print("state ",i,":\n")
    m=secant(a,b,1e-6,s,150,(i+0.5)*2*π)
    b=m
    a=m-0.02
    print(m,"\t value of s(e)= ",s(150,m,0),"\n\n")
    push!(O2,m)
    i=i+1
end
while(i<10)
    print("state ",i,":\n")
    m=secant(-0.34,-0.33,1e-6,s,150,(i+0.5)*2*π)
    print(m,"\t value of s(e)= ",s(150,m,0),"\n\n")
    push!(O2,m)
    i=i+1
end
while(i<13)
    print("state ",i,":\n")
    m=secant(-0.18,-0.15,1e-6,s,150,(i+0.5)*2*π)
    print(m,"\t value of s(e)= ",s(150,m,0),"\n\n")
    push!(O2,m)
    i=i+1
end
while(i<15)
    print("state ",i,":\n")
    m=secant(-0.08,-0.05,1e-6,s,150,(i+0.5)*2*π)
    print(m,"\t value of s(e)= ",s(150,m,0),"\n\n")
    push!(O2,m)
    i=i+1
end
while(i<16)
    print("state ",i,":\n")
    m=secant(-0.03,-0.02,1e-6,s,150,(i+0.5)*2*π)
    print(m,"\t value of s(e)= ",s(150,m,0),"\n\n")
    push!(O2,m)
    i=i+1
end
print("No more states for oxygen could be found using this program!\n")
x=range(0.95,stop=6,step=1/6000)
p=plot(x,LJ.(x,0),"-",label="Lennard Jonnes Potential")
i=0
#plotting the states for O2
for y in O2
    x1=range(r_infind(y),stop=r_outfind(y),step=1/600)
    y1=[y for x in x1]
    p=plot(x1,y1,label=string(round(y,digits=3))*" dimensionless units")
    i=i+1
end
PyPlot.legend()
PyPlot.grid("on")

#plotting the states for H2
x=range(0.95,stop=6,step=1/6000)
p=plot(x,LJ.(x,0),"-",label="Lennard Jonnes Potential")
i=0
for y in H2
    x1=range(r_infind(y),stop=r_outfind(y),step=1/600)
    y1=[y for x in x1]
    p=plot(x1,y1,label=string(round(y,digits=3))*" dimensionless units")
    i=i+1
end
PyPlot.legend()
PyPlot.grid("on")

#plotting the states for H2
x=range(0.95,stop=6,step=1/6000)
p=plot(x,LJ.(x,0),"-",label="Lennard Jonnes Potential")
i=0
for y in H2
    x1=range(r_infind(y),stop=r_outfind(y),step=1/600)
    y1=[y for x in x1]
    p=plot(x1,y1,label=string(round(y,digits=3))*" dimensionless units")
    i=i+1
end
PyPlot.legend()
PyPlot.grid("on")
