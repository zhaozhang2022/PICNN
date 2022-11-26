
import torch
import torchvision
import torch.utils.data
import math
import time

import copy

# aa=torch.tensor([1, 2, 3, 4])
# bb=aa.reshape(2,2)

f=open('samplesperm.txt','r')
data = f.readlines()  # 将txt中所有字符串读入data
permvec = list(map(float, data))
f.close()
for i in range(400):
    permvec[i]=permvec[i]*1e-15


class NODE:
    def __init__(self):
        self.x = 0
        self.y = 0
        self.z = 0

class CELL:
    def __init__(self): # 不加self就变成了对所有类对象同时更改
        self.vertices = [-1, -1, -1, -1, -1, -1, -1, -1]
        self.neighbors = [-1, -1, -1, -1, -1, -1]
        self.dx = 0
        self.dy = 0
        self.dz = 0
        self.volume = 0
        self.xc = 0
        self.yc = 0
        self.zc = 0
        self.porosity = 0
        self.kx = 0
        self.ky = 0
        self.kz = 0
        self.trans = [0, 0, 0, 0, 0, 0]
        self.transw = [0, 0, 0, 0, 0, 0]
        self.markbc = -2
        self.press = 0
        self.Sw = 0
        self.markbc_Sw=0
        self.markwell=-1
        self.mobiw=0
        self.mobio=0
        self.mobit=0

print("build Grid")
ddx=5.0
dxvec=[0]
for i in range(0, 20):
    dxvec.append(5)

dyvec=[0]
for i in range(0, 20):
    dyvec.append(5)
dzvec=[0,5]


nx=len(dxvec)-1
ny=len(dyvec)-1
nz=len(dzvec)-1
nodelist=[]
llz = 0
for k in range(0, nz+1):
    llz = llz + dzvec[k]
    lly=0
    for j in range(0, ny+1):
        lly = lly + dyvec[j]
        llx = 0
        for i in range(0, nx+1):
            llx = llx + dxvec[i]
            node=NODE()
            node.x=llx
            node.y=lly
            node.z=llz
            nodelist.append(node)

# build connectivity and neighbors
celllist=[]

for k in range(0, nz):
    for j in range(0, ny):
        for i in range(0, nx):
            id = k * nx * ny + j * nx + i
            nc=id
            cell = CELL()
            if i>0:
                cell.neighbors[0] = nc - 1
            if i<nx-1:
                cell.neighbors[1] = nc + 1
            if j>0:
                cell.neighbors[2] = nc - nx
            if j<ny-1:
                cell.neighbors[3] = nc + nx
            if k>0:
                cell.neighbors[4] = nc - nx*ny
            if k<nz-1:
                cell.neighbors[5] = nc + nx * ny
            i0 = k * (nx + 1) * (ny + 1) + j * (nx + 1) + i
            i1 = k * (nx + 1) * (ny + 1) + j * (nx + 1) + i + 1
            i2 = k * (nx + 1) * (ny + 1) + (j + 1) * (nx + 1) + i
            i3 = k * (nx + 1) * (ny + 1) + (j + 1) * (nx + 1) + i + 1
            i4 = (k + 1) * (nx + 1) * (ny + 1) + j * (nx + 1) + i
            i5 = (k + 1) * (nx + 1) * (ny + 1) + j * (nx + 1) + i + 1
            i6 = (k + 1) * (nx + 1) * (ny + 1) + (j + 1) * (nx + 1) + i
            i7 = (k + 1) * (nx + 1) * (ny + 1) + (j + 1) * (nx + 1) + i + 1
            cell.dx = nodelist[i1].x - nodelist[i0].x
            cell.dy = nodelist[i2].y - nodelist[i0].y
            cell.dz = nodelist[i4].z - nodelist[i0].z
            cell.vertices[0] = i0
            cell.vertices[1] = i1
            cell.vertices[2] = i2
            cell.vertices[3] = i3
            cell.vertices[4] = i4
            cell.vertices[5] = i5
            cell.vertices[6] = i6
            cell.vertices[7] = i7
            cell.xc = 0.125 * (nodelist[i0].x+nodelist[i1].x+nodelist[i2].x+nodelist[i3].x+nodelist[i4].x+nodelist[i5].x+nodelist[i6].x+nodelist[i7].x)
            cell.yc = 0.125 * (nodelist[i0].y + nodelist[i1].y + nodelist[i2].y + nodelist[i3].y + nodelist[i4].y + nodelist[i5].y + nodelist[i6].y + nodelist[i7].y)
            cell.zc = 0.125 * (nodelist[i0].z + nodelist[i1].z + nodelist[i2].z + nodelist[i3].z + nodelist[i4].z + nodelist[i5].z + nodelist[i6].z + nodelist[i7].z)
            cell.volume=cell.dx*cell.dy*cell.dz
            celllist.append(cell)

cellvolume=celllist[0].volume
ncell=len(celllist)

print("define properties")
mu_o = 1.8e-3
mu_w = 1e-3
chuk = 15e-15
poro = 0.2
Siw=0.2
Bo = 1.2
Bw = 1.0
Cr = 10 * 1e-6 / 6894
Cw = 4 * 1e-6 / 6894
Co = 100 * 1e-6 / 6894
p_init = 20e6
p_e = 20e6



print("set properties to grid and initial conditions")
for i in range(0, ncell):
    celllist[i].porosity=poro
    # celllist[i].kx = permvec[i]
    # celllist[i].ky = permvec[i]
    # celllist[i].kz = permvec[i]
    celllist[i].kx = chuk
    celllist[i].ky = chuk
    celllist[i].kz = chuk
    celllist[i].Sw = Siw
    celllist[i].press=p_init


print("set well conditions")
celllist[0].markwell = 0  # 注水井
celllist[0].markbc = -1
celllist[0].markbc_Sw = 1
celllist[0].Sw = 1
celllist[ncell - 1].markwell = 1  # 生产井
celllist[ncell - 1].markbc = -1
celllist[nx - 1].markbc = 1
celllist[nx - 1].press = p_e
celllist[nx * nx - nx].markbc = 1
celllist[nx * nx - nx].press = p_e



print("mobility function")
def computemobi():
    for ie in range(0, ncell):
        sw=celllist[ie].Sw
        a=(1-sw)/(1-Siw)
        b=(sw-Siw)/(1-Siw)
        kro=a*a*(1-b*b)
        krw=b*b*b*b
        vro=kro/mu_o/Bo
        vrw=krw/mu_w/Bw
        celllist[ie].mobio=vro
        celllist[ie].mobiw=vrw
        celllist[ie].mobit=vro+vrw

print("transmissibility function")
def computetrans():
    for ie in range(0, ncell):
        for j in range(0, 4):
            je = celllist[ie].neighbors[j]
            if je >= 0:
                mt1=celllist[ie].mobit
                mt2=celllist[je].mobit
                k1 = celllist[ie].kx
                k2 = celllist[je].kx
                t1=mt1*k1*ddx*ddx/(ddx/2)
                t2=mt2*k2*ddx*ddx/(ddx/2)
                tt=1 / (1 / t1 + 1 / t2)
                celllist[ie].trans[j] =tt

print("neighbour tensor")

neiborvec_w = torch.zeros(ncell).type(torch.long).cuda()
neiborvec_e = torch.zeros(ncell).type(torch.long).cuda()
neiborvec_s = torch.zeros(ncell).type(torch.long).cuda()
neiborvec_n = torch.zeros(ncell).type(torch.long).cuda()


for ie in range(ncell):
    neibor_w = celllist[ie].neighbors[0]
    neibor_e = celllist[ie].neighbors[1]
    neibor_s = celllist[ie].neighbors[2]
    neibor_n = celllist[ie].neighbors[3]
    if neibor_w < 0:
        neibor_w = 0
    if neibor_e < 0:
        neibor_e = 0
    if neibor_s < 0:
        neibor_s = 0
    if neibor_n < 0:
        neibor_n = 0
    neiborvec_w[ie] = neibor_w
    neiborvec_e[ie] = neibor_e
    neiborvec_n[ie] = neibor_n
    neiborvec_s[ie] = neibor_s

print("trans tensor")
transvec_n = torch.zeros(ncell).cuda()
transvec_s = torch.zeros(ncell).cuda()
transvec_e = torch.zeros(ncell).cuda()
transvec_w = torch.zeros(ncell).cuda()

def transcell2tensor():
    for ie in range(ncell):
        transvec_w[ie] = celllist[ie].trans[0]
        transvec_e[ie] = celllist[ie].trans[1]
        transvec_s[ie] = celllist[ie].trans[2]
        transvec_n[ie] = celllist[ie].trans[3]

output_size = ncell

print("define NN model, criterion, optimizer and scheduler")

class CNN(torch.nn.Module):
    def __init__(self, output_size):
        super(CNN, self).__init__()
        self.conv1 = torch.nn.Sequential(
            torch.nn.Conv2d(1, 25, kernel_size=(3,3), stride=(1,1), padding=1), # 20*20 -> 25*20*20
            torch.nn.BatchNorm2d(25),
            torch.nn.ReLU(),
            # torch.nn.Dropout(),
            torch.nn.MaxPool2d(kernel_size=2, stride=2)  # # ->25*10*10
        )

        self.conv2 = torch.nn.Sequential(
            torch.nn.Conv2d(25, 50, kernel_size=(3,3), stride=(1,1), padding=1), # 25*10*10 -> 50*10*10
            torch.nn.BatchNorm2d(50),
            torch.nn.ReLU(),
            # torch.nn.Dropout(),
            torch.nn.MaxPool2d(kernel_size=2, stride=2)  # ->50*5*5
        )

        self.fc = torch.nn.Sequential(
            torch.nn.Linear(50 * 5 * 5, output_size),
            # torch.nn.ReLU(),
            # torch.nn.Linear(600, output_size)
        )

    def forward(self, x):
        x = self.conv1(x)
        x = self.conv2(x)
        x = x.reshape(x.size(0), -1)   # x: tensor (time steps (100), 50*3*3)
        x = self.fc(x) # x: tensor (time steps (100), output size)
        return x


# 实例化模型
model = CNN(output_size)
model.cuda()
# 设置损失函数和优化器
learning_rate = 0.01
criterion = torch.nn.MSELoss()
# criterion = torch.nn.L1Loss()
# criterion = torch.nn.L1Loss(reduction='none')
optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)
#optimizer = torch.optim.SGD(model.parameters(), lr=learning_rate, momentum=0.99)
num_epochs = 2000
floss = open('alllosslowest_cnn_homo.txt','w')
writepwf=open('pwfproducer_cnn_homo.txt','w')
writepe=open('peproducer_cnn_homo.txt','w')
f = open('result_CNN_homo.vtk','w')
ftime=open('trainingtime_cnn_homo20x20.txt','w')
print("训练模型")

def pdeimplicit(p, presslast, alphavec, qtvec):  #bound指的是确保effectively输出是在0到1之间起作用, 因此用abs。
    pp = torch.zeros_like(p).cuda()
    pp[:]=abs(p[:])*p_init
    pp[nx-1]=p_e
    pp[nx*nx-nx]=p_e
    pde1=torch.zeros_like(p).cuda()
    pde1[:] = pde1[:] - transvec_w[:] * (pp[neiborvec_w[:]] - pp[:])
    pde1[:] = pde1[:] - transvec_e[:] * (pp[neiborvec_e[:]] - pp[:])
    pde1[:] = pde1[:] - transvec_s[:] * (pp[neiborvec_s[:]] - pp[:])
    pde1[:] = pde1[:] - transvec_n[:] * (pp[neiborvec_n[:]] - pp[:])
    pde1[:] = pde1[:] - qtvec[:] + (pp[:]-presslast[:])*alphavec[:]
    return pde1

print("construct input tensor")
inputtensor = torch.ones(1, 1, ny, nx).cuda()
resultout=torch.ones(ncell).cuda()
presslast=torch.ones(ncell).cuda()
qtvec=torch.zeros(ncell).cuda()
for i in range(ncell):
    resultout[i]=(celllist[i].press)/p_init
    presslast[i]=celllist[i].press
qin=5.0/86400
qtvec[0]= qin
qtvec[nx*nx-1]=-qin
print("Time Iteration")
nt = 500
dt=36000
alphavec = torch.zeros(ncell).cuda()
re = 0.14*(ddx*ddx + ddx*ddx)**0.5
SS=3
rw=0.05
writepe.write("%e\n" % p_init)
writepwf.write("%e\n" % p_init)

totaltime=0
for t in range(nt):
    print('Time is ', t)
    if t>0:
        num_epochs=500
    lowestloss = 10000000000
    computemobi()
    computetrans()
    transcell2tensor()
    for ie in range(ncell):
        alphavec[ie]=celllist[ie].porosity*(1-celllist[ie].Sw)*(Cr+Co)/Bo+celllist[ie].porosity*celllist[ie].Sw*(Cr+Cw)/Bw
        alphavec[ie] = alphavec[ie]*celllist[ie].volume/dt
    inputtensor[0][0]=resultout.reshape(ny, nx)
    print("CNN Implicit Solver")
    starttime = time.time()
    for epoch in range(num_epochs):
        outputtensor = model(inputtensor)  # on gpu
        resultnext = outputtensor[0].clone().detach()
        diff=pdeimplicit(outputtensor[0], presslast, alphavec, qtvec)
        # 计算损失并利用反向传播计算损失对各参数梯度
        loss = criterion(diff, diff * 0)
        optimizer.zero_grad()
        loss.backward()
        # loss.backward(torch.ones_like(loss))
        optimizer.step()
        # scheduler.step()
        if loss < lowestloss:
            lowestloss = loss
            resultout = resultnext
        if epoch % 100 == 0:
            print('epoch is ', epoch, 'lowestloss is: ', lowestloss, '\n')
            floss.write("%e\n" % lowestloss)
    endtime = time.time()
    simtime = endtime - starttime
    totaltime=totaltime+simtime
    for ie in range(ncell):
        resultout[ie]=abs(resultout[ie])
    for ie in range(ncell):
        celllist[ie].press = resultout[ie] * p_init
        celllist[nx-1].press = p_e
        celllist[nx*nx-nx].press = p_e
    writepe.write("%e\n" % celllist[nx*nx-1].press)
    PI = 2 * 3.14 * ddx * celllist[nx * nx - 1].kx * celllist[nx * nx - 1].mobit / (math.log(re / rw) + SS)
    pwf=celllist[nx*nx-1].press-qin/PI
    writepwf.write("%e\n" % pwf)
    #update saturation
    for ie in range(ncell):
        if celllist[ie].markbc_Sw==0:
            tfluxsw=0
            tfluxin=0
            pi=celllist[ie].press
            for i in range(4):
                je=celllist[ie].neighbors[i]
                if je>=0:
                    pj=celllist[je].press
                    if pj>pi:
                        fluxin=(pj-pi)*celllist[ie].trans[i]
                        tfluxin += fluxin
                        tfluxsw += fluxin*celllist[je].mobiw/celllist[je].mobit
            tfluxout=-tfluxin
            tfluxsw += tfluxout*celllist[ie].mobiw/celllist[ie].mobit
            if ie==0:
                tfluxsw += qin
            sw=celllist[ie].Sw
            tfluxsw += -(pi-presslast[ie])/dt*poro*sw/Bw*(Cr+Cw)
            celllist[ie].Sw = celllist[ie].Sw + tfluxsw*dt/poro/celllist[ie].volume*Bw
    for ie in range(ncell):
        presslast[ie] = celllist[ie].press


ftime.write("%0.3f\n" % totaltime)
print(totaltime)
print("output to vtk")



f.write("# vtk DataFile Version 2.0\n")
f.write( "Unstructured Grid\n")
f.write( "ASCII\n")
f.write("DATASET UNSTRUCTURED_GRID\n")
f.write("POINTS %d double\n" % (len(nodelist)))
for i in range(0, len(nodelist)):
    f.write("%0.3f %0.3f %0.3f\n" % (nodelist[i].x, nodelist[i].y, nodelist[i].z))
f.write("\n")
f.write("CELLS %d %d\n" % (len(celllist), len(celllist)*9))
for i in range(0, len(celllist)):
    f.write("%d %d %d %d %d %d %d %d %d\n" % (8, celllist[i].vertices[0], celllist[i].vertices[1], celllist[i].vertices[3], celllist[i].vertices[2], celllist[i].vertices[4], celllist[i].vertices[5], celllist[i].vertices[7], celllist[i].vertices[6]))
f.write("\n")
f.write("CELL_TYPES %d\n" % (len(celllist)))
for i in range(0, len(celllist)):
    f.write("12\n")
f.write("\n")
f.write("CELL_DATA %d\n" % (len(celllist)))
f.write("SCALARS Pressure double\n")
f.write("LOOKUP_TABLE default\n")
for i in range(0, len(celllist)):
    f.write("%0.3f\n" % (celllist[i].press/10**6))
f.write("\n")
f.write("SCALARS Saturation double\n")
f.write("LOOKUP_TABLE default\n")
for i in range(0, len(celllist)):
    f.write("%0.3f\n" % (celllist[i].Sw))
f.close()







