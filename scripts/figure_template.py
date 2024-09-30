import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from IPython.display import HTML, display, FileLink
from base64 import b64encode, b64decode
from io import StringIO, BytesIO
from wand.image import Image as WImage
from wand.display import display as wdisplay
from matplotlib.ticker import FixedFormatter, AutoLocator, AutoMinorLocator, MaxNLocator
from matplotlib.backends.backend_pdf import PdfPages, PdfFile
from contextlib import contextmanager
import string
import uuid

__all__ = [
    'cm_to_inch', 'inch_to_cm', 'image_link', 'auto_xticklabels', 'auto_yticklabels',
    'legend', 'plot_demo', 'download_button', 'download_figure', 
    'embed_pdf_figure', 'embed_pdf_pages',  'embed_pdf_data',
    'display_dataframe', 'setup_theme', 'std_plot'
]

def cm_to_inch(width, height):
    return width/2.54, height/2.54

def inch_to_cm(width, height):
    return width*2.54, height*2.54

def image_link(filename):
    plt.savefig(filename)
    display(WImage(filename=filename))
    plt.close()
    #display(FileLink(filename))
    display(HTML('<a href="{0}" download="{0}" target="_blank">{0}</a>'.format(filename)))

def auto_xticklabels(ax, fontdict=None, nbins='auto', steps=None, **kwargs):
    if fontdict is None:
        fontdict={'fontweight': 'bold'}
    if steps is None:
        steps = [1, 2, 2.5, 5, 10]
    locator = MaxNLocator(nbins=nbins, steps=steps, **kwargs)
    locator.set_axis(ax.xaxis)
    ax.set_xticks(locator())
    ax.set_xticklabels(locator(), fontdict=fontdict)
    
def auto_yticklabels(ax, fontdict=None, nbins='auto', steps=None, **kwargs):
    if fontdict is None:
        fontdict={'fontweight': 'bold'}
    if steps is None:
        steps = [1, 2, 2.5, 5, 10]
    locator = MaxNLocator(nbins=nbins, steps=steps, **kwargs)
    locator.set_axis(ax.yaxis)
    ax.set_yticks(locator())
    ax.set_yticklabels(locator(), fontdict=fontdict)
    
def legend(ax, **kwargs):
    l = ax.legend(**kwargs)
    l.get_title().set_fontweight('bold')
    l.get_title().set_fontsize(12)
    return l

def plot_demo():
    x = np.linspace(0, 10, 1000)
    fig, ax = plt.subplots(figsize=cm_to_inch(12, 9))
    ax.plot(x, np.sin(np.pi*x), label='sin(x)', linewidth=1.5)
    ax.plot(x, np.cos(np.pi*x), label='cos(y)', linewidth=1.5)
    ax.set_xlabel('xlabel')
    ax.set_ylabel('ylabel')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    ax.spines['left'].set_position(('outward', 10))
    ax.spines['bottom'].set_position(('outward', 10))
    auto_xticklabels(ax)
    auto_yticklabels(ax)
    ax.set_ylim(-1, 1)
    ax.set_xlim(0, 10)
    
    ax.set_title('Font size demo')
    legend(ax, title='Curve', bbox_to_anchor=(1.04,0.5), loc="center left", borderaxespad=0)
    fig.tight_layout()

def download_button(data, filename=None, mime_type='application/octet-stream', label='Download', encoding='utf-8'):
    if filename is None:
        filename = 'download'
    data = 'data:{mime_type};base64,'.format(mime_type=mime_type) \
        + str(b64encode(data), encoding=encoding)
    button_id = 'button_{}'.format(np.random.randint(1000000000))
    display(HTML(r'<input type="button" id="{0}" value="Download">'.format(button_id)))
    display(HTML('''<script>
    document.getElementById("{button_id}").addEventListener("click", function(event){{
        var filename = "{filename}";
        var data = "{data}";
        const element = document.createElement('a');
        element.setAttribute('href', data);
        element.setAttribute('download', filename);
        element.style.display = 'none';
        document.body.appendChild(element);
        element.click();
        document.body.removeChild(element);
    }});
</script>'''.format(button_id=button_id, filename=filename, data=data)))

def download_figure(filename, format=None):
    if format is None:
        format = filename.split('.')[-1]
    data = BytesIO()
    plt.savefig(data, format=format)
    data = data.getvalue()
    download_button(data, filename=filename)

"""
def embed_pdf_figure(width=640, height=480, title='Image'):
    data = BytesIO()
    plt.savefig(data, format='pdf', metadata={'Title': title})
    data = data.getvalue()
    data = 'data:application/pdf;base64,'+ str(b64encode(data), encoding='utf-8')
    display(HTML('<object width="{}" height="{}" data="{}"></object>'.format(width, height, data)))
    plt.close()
"""

def render_pdf_html(width, height, title, data):
    return HTML(string.Template('''
<div id="container-${id}"></div>
<script>
{
    let data = "${data}";
    let div = document.getElementById("container-${id}");
    let o = document.createElement("object");
    o.setAttribute("width", "${width}");
    o.setAttribute("height", "${height}");
    o.setAttribute("data", data);
    div.appendChild(o);
    div.appendChild(document.createElement("br"));
    let button = document.createElement("input");
    button.setAttribute("type", "button");
    button.setAttribute("id", "button-${id}");
    button.setAttribute("value", "Download");
    button.addEventListener("click", function(event){
        let filename = "${title}.pdf";
        const element = document.createElement('a');
        element.setAttribute("href", data);
        element.setAttribute("download", filename);
        element.style.display = "none";
        document.body.appendChild(element);
        element.click();
        document.body.removeChild(element);
    });
    div.appendChild(button);
}
</script>''').substitute(id=uuid.uuid1().hex, width=width, height=height, title=title, data=data))


def embed_pdf_figure(width=960, height=480, title='Image'):
    data = BytesIO()
    plt.savefig(data, format='pdf', metadata={'Title': title})
    data = data.getvalue()
    data = 'data:application/pdf;base64,'+ str(b64encode(data), encoding='utf-8')
    html = render_pdf_html(width=width, height=height, title=title, data=data)
    display(html)
    plt.close()

@contextmanager
def embed_pdf_data(width=640, height=480, title='Image'):
    try:
        data = BytesIO()
        yield data
    finally:
        data = data.getvalue()
        data = 'data:application/pdf;base64,'+ str(b64encode(data), encoding='utf-8')
        display(render_pdf_html(width=width, height=height, data=data, title=title))
        plt.close()

def embed_pdf_grid(g, width=960, height=480, title='Image'):
    data = BytesIO()
    g.savefig(data, format='pdf', metadata={'Title': title})
    data = data.getvalue()
    data = 'data:application/pdf;base64,'+ str(b64encode(data), encoding='utf-8')
    html = render_pdf_html(width=width, height=height, title=title, data=data)
    display(html)
    plt.close()

@contextmanager
def embed_pdf_pages(width=960, height=480, title='Image'):
    '''Embed PDF with multiple pages in Jupyter
    Example:
        with embed_pdf_pages() as pdf:
            for T in range(1, 5):
                fig, ax = plt.subplots()
                x = np.linspace(0, 10)
                y = np.sin(2*np.pi*x)
                ax.plot(x, y)
                pdf.savefig(fig)
                plt.close()
    '''
    data = BytesIO()
    try:
        pdf = PdfPages(data, metadata={'Title': title})
        yield pdf
    finally:
        pdf.close()
        data = data.getvalue()
        data = 'data:application/pdf;base64,'+ str(b64encode(data), encoding='utf-8')
        display(render_pdf_html(width=width, height=height, data=data, title=title))
        plt.close()

def display_dataframe(df, filename=None, encoding='utf-8', format='csv', type='button', **kwargs):
    display(df)
    if filename is None:
        filename = "dataframe"
    if format == 'csv':
        data = bytes(df.to_csv(**kwargs), encoding=encoding)
        mime_type = 'text/csv'
        filename = filename + '.csv'
    elif format == 'tsv':
        data = bytes(df.to_csv(sep='\t', **kwargs), encoding=encoding)
        mime_type = 'text/plain'
        filename = filename + '.txt'
    elif format == 'excel':
        data = BytesIO()
        df.to_excel(pd.ExcelWriter(data))
        data = data.getvalue()
        filename = filename + '.xlsx'
        mime_type = 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet'
    else:
        raise ValueError('unknown file format: {}'.format(format))
    data = 'data:{mime_type};base64,'.format(mime_type=mime_type) + str(b64encode(data), encoding=encoding)
    if type == 'hyperlink':
        display(HTML('<a href="{data}" download={filename} target="_blank">{filename}</a>'.format(
            mime_type=mime_type, filename=filename, data=data)))
    elif type == 'button':
        button_id = 'button_{}'.format(uuid.uuid1())
        display(HTML(r'<input type="button" id="{0}" value="Download">'.format(button_id)))
        display(HTML('''<script>
    document.getElementById("{button_id}").addEventListener("click", function(event){{
        var filename = "{filename}";
        var data = "{data}";
        const element = document.createElement('a');
        element.setAttribute('href', data);
        element.setAttribute('download', filename);
        element.style.display = 'none';
        document.body.appendChild(element);
        element.click();
        document.body.removeChild(element);
    }});
</script>'''.format(button_id=button_id, filename=filename, data=data)))

def setup_theme():
    sns.set_style('white')
    #plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['figure.dpi'] = 120
    plt.rcParams['font.size'] = 10
    plt.rcParams['axes.titlesize'] = 14
    plt.rcParams['axes.labelsize'] = 12
    plt.rcParams['axes.labelweight'] = 'bold'
    plt.rcParams['xtick.labelsize'] = 12
    plt.rcParams['ytick.labelsize'] = 12
    plt.rcParams['xtick.major.size'] = 5
    plt.rcParams['xtick.major.width'] = 1
    plt.rcParams['ytick.major.size'] = 5
    plt.rcParams['ytick.major.width'] = 1
    plt.rcParams['axes.titleweight'] = 'bold'
    plt.rcParams['legend.fontsize'] = 12
    plt.rcParams['xtick.bottom'] = True
    plt.rcParams['ytick.left'] = True
    plt.rcParams['axes.linewidth'] = 1.0
    plt.rcParams['patch.linewidth'] = 1.0

def std_plot(ax,xlabel,ylabel,title=None,
             legendtitle=None,bbox_to_anchor=None,
             labelspacing=1.2,borderpad=0,handletextpad=0,legendsort=True,markerscale=None,
             xlim=None,ylim=None,
             xbins=None,ybins=None,
             cbar=None,cbarlabel=None,
             moveyaxis=False,sns=False,left=True):
    fonttitle = {'family':'Arial',
                  'weight' : 'normal', 
                  'size' : 8}
    fontlabel = {'family':'Arial',
                    'weight' : 'normal', 
                    'size' : 6.5}
    fontticklabel = {'family':'Arial',
                    'weight' : 'normal', 
                    'size' : 6.5}
    fontlegend = {'family':'Arial',
                    'weight' : 'normal', 
                #'linewidth':0.5,
                    'size' : 6.5}
    fontcbarlabel = {'family':'Arial',
                    'weight' : 'normal', 
                    #'Rotation' : 270,
                    #'labelpad' : 25,
                    'size' : 6.5}
    fontcbarticklabel = {'family':'Arial',#Helvetica
                    'weight' : 'normal', 
                    'size' : 5.5}
    #plt.yticks([0,0.2,0.4,0.6,0.8,1.0],['0','0.2','0.4','0.6','0.8','1.0'])
    plt.draw()
    #plt.figure(linewidth=30.5)
    if xlim is not None:  
        ax.set(xlim=xlim)
    if ylim is not None:
        ax.set(ylim=ylim)
    if xbins is not None:
        locator = MaxNLocator(nbins=xbins)
        locator.set_axis(ax.xaxis)
        ax.set_xticks(locator())
    if ybins is not None:
        locator = MaxNLocator(nbins=ybins)
        locator.set_axis(ax.yaxis)
        ax.set_yticks(locator())
    plt.draw()
    ax.set_xticks(ax.get_xticks())
    ax.set_yticks(ax.get_yticks())
    ax.set_xlabel(xlabel,fontdict = fontlabel,labelpad=5.5)
    ax.set_ylabel(ylabel,fontdict = fontlabel,labelpad=5.5)
    ax.set_xticklabels(ax.get_xticklabels(),fontticklabel)
    ax.set_yticklabels(ax.get_yticklabels(),fontticklabel)

    if moveyaxis is True:
        #fontticklabel 
        ax.spines['left'].set_position(('data',0))
    ax.spines['left'].set_visible(left)
    ax.spines['right'].set_visible(not left)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_linewidth(0.5)
    ax.spines['bottom'].set_linewidth(0.5)
    ax.spines['left'].set_linewidth(0.5)
    ax.spines['bottom'].set_color('k')
    ax.spines['left'].set_color('k')
    ax.spines['right'].set_color('k')
    
    ax.tick_params(direction='out', pad=2)
    #ax.spines['bottom']._edgecolor="#000000"
    #ax.spines['left']._edgecolor="#000000"
    if title is not None:
        ax.set_title(title,fontdict = fonttitle)
    if legendtitle is not None:
        #if legendloc is None:
        #    legendloc="best"
        legend = ax.legend(title=legendtitle,prop=fontlegend,
                      bbox_to_anchor=bbox_to_anchor,
                      labelspacing=labelspacing,borderpad=borderpad,handletextpad=handletextpad,
                      edgecolor="#000000",fancybox=False,markerscale=markerscale)
        ax.legend_.get_frame()._linewidth=0.5
        legend.get_title().set_fontweight('normal')
        legend.get_title().set_fontsize(6.5)
        if legendsort is True:
            # h: handle l:label
            h,l = ax.get_legend_handles_labels()
            l,h = zip(*sorted(zip(l,h), key=lambda t: int(t[0]))) 
            legend = ax.legend(h,l,title=legendtitle,prop=fontlegend,
                      bbox_to_anchor=bbox_to_anchor,
                      labelspacing=labelspacing,borderpad=borderpad,handletextpad=handletextpad,
                      edgecolor="#000000",fancybox=False,markerscale=markerscale)
            ax.legend_.get_frame()._linewidth=0.5
            legend.get_title().set_fontweight('normal')
            legend.get_title().set_fontsize(6.5)
        if sns is True:
            h,l = ax.get_legend_handles_labels()
            #l,h = zip(*sorted(zip(l,h), key=lambda t: int(t[0]))) 
            legend = ax.legend(h[1:],l[1:],title=legendtitle,prop=fontlegend,
                      bbox_to_anchor=bbox_to_anchor,
                      labelspacing=labelspacing,borderpad=borderpad,handletextpad=handletextpad,
                      edgecolor="#000000",fancybox=False,markerscale=markerscale)
            ax.legend_.get_frame()._linewidth=0.5
            legend.get_title().set_fontweight('normal')
            legend.get_title().set_fontsize(6.5)

    if cbar is not None:
        #locator, formatter = cbar._get_ticker_locator_formatter()
        #ticks, ticklabels, offset_string = cbar._ticker(locator, formatter)
        #cbar.ax.spines['top'].set_visible(False)
        #cbar.ax.spines['right'].set_visible(False)
        #cbar.ax.spines['bottom'].set_visible(False)
        #cbar.ax.spines['left'].set_visible(False)
        cbar.ax.tick_params(direction='out', pad=3,width=0,length=0)
        cbar.set_label(cbarlabel,fontdict = fontcbarlabel,Rotation=270,labelpad=7.5)
        cbar.ax.set_yticks(cbar.ax.get_yticks())
        cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(),fontcbarticklabel)
    return ax