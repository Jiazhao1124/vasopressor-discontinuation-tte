import matplotlib as mpl
mpl.use('Agg')
from matplotlib import font_manager as fm
import os
for f in ['Regular','Bold']:
    _p = f'/usr/share/fonts/truetype/liberation/LiberationSans-{f}.ttf'
    if os.path.exists(_p):
        fm.fontManager.addfont(_p)  # metric-compatible Arial substitute on Linux
mpl.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','Liberation Sans'],'font.size':8,'axes.linewidth':0.6,
 'xtick.major.width':0.6,'ytick.major.width':0.6,'xtick.major.size':3,'ytick.major.size':3,'xtick.direction':'out','ytick.direction':'out',
 'svg.fonttype':'none','pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
BLUE='#0072B5'; RED='#BC3C29'; GREY='#4D4D4D'; LIGHT='#9E9E9E'
def save(fig,name):
    for ext,kw in [('png',{'dpi':600}),('tiff',{'dpi':600,'pil_kwargs':{'compression':'tiff_lzw'}}),('pdf',{}),('svg',{})]:
        fig.savefig(f'{name}.{ext}',**kw)
