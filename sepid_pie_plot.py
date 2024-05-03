from sepid import *

''' designed for (nested) pie plot '''

plt.close('all')
datadir = '/scratch/sepid/DATA/AR/plage/2016.09.15/'
savedir = datadir+'OUTPUT/'
outdir = '/home/seki2695/OUTPUT/inv/'

# data counts
a = 24
b = 60
c = 60
d = 27
cmap = plt.get_cmap("tab20c")

# outer pie
rad_out = 1.
labels_out = ['similar\nstructure', 'no clear\ncorrespondence']
sizes_out = [a+b , c+d]
colors_out =['tomato', 'gray']
exp_out = [0.01,0.01]

# inner pie
rad_in = 0.66
labels_in = ['(a)','(b)','(c)','(d)']
sizes_in = [a,b,c,d]
colors_in = ['coral', 'lightsalmon', 'lightgrey','darkgrey']
exp_in = [0.01,0.01,0.01,0.01]

# center circle
rad_c = 0.33
 
# outer pie Plot
patches, texts, autotexts = plt.pie(sizes_out, labels=labels_out, radius = rad_out, colors=colors_out, startangle=90,autopct='%1.1f%%',pctdistance = rad_in+(rad_out-rad_in)/2, labeldistance = rad_out+0.1,frame=True, textprops = {'fontsize':12})#,explode = exp_out, shadow = True)
texts[0].set_fontsize(11)
texts[1].set_fontsize(11)
# inner pie plot
plt.pie(sizes_in,labels = labels_in,colors=colors_in,radius=rad_in,startangle=90,autopct='%1.1f%%',labeldistance = rad_c,pctdistance = rad_in+0.1,textprops = {'fontsize':10.5})#,explode = exp_in, shadow = True)

# inner circle plot
centre_circle = plt.Circle((0,0),rad_c,color='black', fc='white',linewidth=0)
fig = plt.gcf()
fig.set_size_inches(5.75,3.5)
fig.gca().add_artist(centre_circle)
 
plt.axis('equal')

plt.subplots_adjust(left = 0.17,
                    bottom = 0.0,
                    right = 0.75,
                    top = 1.0,
                    wspace = 0.11,
                    hspace = 0.11
)
plt.show()
outname = outdir+'pie_plot_new.pdf'
plt.savefig(outname, quality = 100)
print 'file saved to '+ outname
