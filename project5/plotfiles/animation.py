import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
plt.rcParams['animation.ffmpeg_path'] = '/usr/bin/ffmpeg'
import pyarma as pa
import sys

#
# Let's generate a dummy time series for a function z(x,y,t)
#
filename = sys.argv[1]
number_of_snaps = int(sys.argv[2])

snaps = []
for i in range(number_of_snaps):
    snaps.append(float(sys.argv[i+3]))

# Set up a 2D xy grid
h = 0.005
x_points = np.arange(0, 1+h, h)
y_points = np.arange(0, 1+h, h)
x, y = np.meshgrid(x_points, y_points, sparse=True)

# Array of time points
dt = 0.005
t_points = np.arange(0, 1+dt, dt)

A = pa.cube() #Create pa.mat object (just as arma::mat in C++)
A.load("./datafiles/"+str(filename)) #Load the content of the matrix you saved into your Python program.
# A function for a Gaussian that is travelling
# in the x direction and broadening as time passes


# Fill z_data_list with f(x,y,t)
z_data_list = []
snapshot_index_list = []
c = 0
for t in t_points:
    z_data = np.rot90(np.array(A[pa.single_slice, c]))
    c += 1
    z_data_list.append(z_data)

    #Finds if the timestap matches a snapshot an appends index
    if any(t== t_snap for t_snap in snaps):
        snapshot_index_list.append(int(c-1))





#
# Now the list z_data_list contains a series of "frames" of z(x,y,t),
# where each frame can be plotted as a 2D image using imshow. Let's
# animate it!
#

# Some settings
fontsize = 12
t_min = t_points[0]
x_min, x_max = x_points[0], x_points[-1]
y_min, y_max = y_points[0], y_points[-1]

# Create figure
fig = plt.figure()
ax = plt.gca()




#plt.savefig('./figures'+filename+'_firstframe.pdf')

# Plots for problem7.1, will not be used to solve problem.
#norm = matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(z_data_list[22]))
#img7_1_1 = ax.imshow(z_data_list[22], extent=[x_min,x_max,y_min,y_max], cmap=plt.get_cmap("viridis"), norm=norm)
#plt.savefig('./figures'+filename+'_time0_11.pdf')
#norm = matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(z_data_list[44]))
#img7_1_1 = ax.imshow(z_data_list[44], extent=[x_min,x_max,y_min,y_max], cmap=plt.get_cmap("viridis"), norm=norm)
#plt.savefig('./figures'+filename+'_time0_22.pdf')
#norm = matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(z_data_list[140]))
#img7_1_2 = ax.imshow(z_data_list[140], extent=[x_min,x_max,y_min,y_max], cmap=plt.get_cmap("viridis"), norm=norm)
#plt.savefig('./figures'+filename+'_time0_70.pdf')
#norm = matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(z_data_list[180]))
#img7_1_3 = ax.imshow(z_data_list[180], extent=[x_min,x_max,y_min,y_max], cmap=plt.get_cmap("viridis"), norm=norm)
#plt.savefig('./figures'+filename+'_time0_90.pdf')
#norm = matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(z_data_list[200]))
#img7_1_4 = ax.imshow(z_data_list[200], extent=[x_min,x_max,y_min,y_max], cmap=plt.get_cmap("viridis"), norm=norm)
#plt.savefig('./figures'+filename+'_time1_0.pdf')


#-----Setting up plot format------
norm = matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(z_data_list[0]))
img = ax.imshow(z_data_list[0], extent=[x_min,x_max,y_min,y_max], cmap=plt.get_cmap("viridis"), norm=norm)


# Axis labels
plt.xlabel("x", fontsize=fontsize)
plt.ylabel("y", fontsize=fontsize)
plt.xticks(fontsize=fontsize)
plt.yticks(fontsize=fontsize)

# Add a colourbar
cbar = fig.colorbar(img, ax=ax)
cbar.set_label("z(x,y,t)", fontsize=fontsize)
cbar.ax.tick_params(labelsize=fontsize)

# Add a text element showing the time
time_txt = plt.text(0.95, 0.95, "t = {:.3e}".format(t_min), color="white",
                    horizontalalignment="right", verticalalignment="top", fontsize=fontsize)
#-----Setting up plot format(end)------


filename = filename.split('.')[0]
#Problem 8.1

for ind in snapshot_index_list:
    norm = matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(z_data_list[ind]))
    img = ax.imshow(z_data_list[ind], extent=[x_min,x_max,y_min,y_max], cmap=plt.get_cmap("viridis"), norm=norm)
    img.set_norm(norm)

    time_txt.set_text("t = {:.3e}".format(t_points[ind]))
    plt.savefig('./figures/'+filename+'_time'+str(t_points[ind])+'.pdf')




img = ax.imshow(z_data_list[0], extent=[x_min,x_max,y_min,y_max], cmap=plt.get_cmap("viridis"), norm=norm)



# Function that takes care of updating the z data and other things for each frame
def animation(i):
    # Normalize the colour scale to the current frame?
    norm = matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(z_data_list[i]))
    img.set_norm(norm)

    # Update z data
    img.set_data(z_data_list[i])

    # Update the time label
    current_time = t_min + i * dt
    time_txt.set_text("t = {:.3e}".format(current_time))

    return img


# Use matplotlib.animation.FuncAnimation to put it all together
anim = FuncAnimation(fig, animation, interval=1, frames=np.arange(0, len(z_data_list), 2), repeat=False, blit=0)

# Run the animation!
plt.show()

# # Save the animation
anim.save('./figures/'+filename+'_animation.gif', writer="ffmpeg", fps=30)
