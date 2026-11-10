"""Final-size Arial 12 pt field plates; native ParaView vector fields retained.
v3 (1 Nov 2026): panel letters only; descriptions moved to captions."""
from pathlib import Path
import json,sys
sys.path.insert(0,str(Path(__file__).resolve().parent))
import fitz
R=Path(__file__).resolve().parents[1];O=R/'figures';P=O/'panels'
M={v['panel']:v for v in json.loads((P/'render_manifest.json').read_text())}
from figure_style import FONT as _F
FONT=str(_F);W=17/2.54*72
COLORS=json.loads((P/'turbo_colors.json').read_text())
def text(p,x,y,s):p.insert_text((x,y),s,fontname='Arial',fontfile=FONT,fontsize=12)
def field(p,stem,rect):
 path=P/(stem+'_vector.svg');assert '<image' not in path.read_text()
 d=fitz.open(path);q=fitz.open('pdf',d.convert_to_pdf())
 p.show_pdf_page(fitz.Rect(rect),q,0,keep_proportion=False)
def chart(p,stem,x,y):
 d=fitz.open(O/(stem+'.pdf'));box=d[0].rect
 # Physical 1:1 placement preserves the chart's 12 pt lettering.
 p.show_pdf_page(fitz.Rect(x,y,x+box.width,y+box.height),d,0)
def scale(p,stem,rect,length,label,xfrac=.06,yfrac=.95,below=False):
 x0,y0,x1,y1=rect;win=M[stem]['world_window'];w=length/(win[1]-win[0])*(x1-x0)
 x=x0+xfrac*(x1-x0);y=y0+yfrac*(y1-y0)
 p.draw_line((x,y),(x+w,y),color=(0,0,0),width=1.5);text(p,x,y+17 if below else y-5,label)
def bar(p,x,y,w,labels,title):
 for i,c in enumerate(COLORS):
  p.draw_rect(fitz.Rect(x+i*w/256,y,x+(i+1)*w/256+.1,y+8),color=None,fill=c)
 p.draw_rect(fitz.Rect(x,y,x+w,y+8),color=(0,0,0),width=.5)
 for f,s in labels:
  width=fitz.Font(fontfile=FONT).text_length(s,fontsize=12)
  text(p,max(3,min(W-width-3,x+f*w-width/2)),y+23,s)
 text(p,x,y-6,title)
def finish(d,name):
 d.save(O/(name+'.pdf'),garbage=4,deflate=True)
 d[0].get_pixmap(matrix=fitz.Matrix(2.8,2.8)).save(O/(name+'.png'))
 (O/(name+'.svg')).write_text(d[0].get_svg_image(text_as_path=False),encoding='utf-8');d.close()

d=fitz.open();p=d.new_page(width=W,height=483)
for stem,y,title in [('G3_speed',25,'(a)'),('G3_pressure',160,'(b)'),('G1_speed',295,'(c)')]:
 text(p,8,y-12,title);rect=(8,y,345,y+94.6);field(p,stem,rect)
 scale(p,stem,rect,20,'20 µm')
 if stem=='G3_speed':bar(p,365,y+30,104,[(0,'0'),(.5,'2'),(1,'4')],'Speed / U')
 if stem=='G3_pressure':
  text(p,361,y+5,'Normalized');text(p,361,y+20,'pressure')
  bar(p,365,y+45,104,[(0,'0'),(.5,'0.5'),(1,'1')],'')
 if stem=='G1_speed':
  win=M[stem]['world_window'];x=lambda v:rect[0]+(v-win[0])/(win[1]-win[0])*(rect[2]-rect[0]);yy=lambda v:rect[3]-(v-win[2])/(win[3]-win[2])*(rect[3]-rect[1])
  box=fitz.Rect(x(37),yy(-3.5),x(44.5),yy(-11.5));p.draw_rect(box,color=(.737,.337,.094),width=.9)
  p.draw_line(box.top_right,(361,311),color=(.737,.337,.094),width=.7)
text(p,359,283,'(d)')
rect=(359,298,472,418.5);field(p,'G1_pocket',rect);scale(p,'G1_pocket',rect,2,'2 µm',.06,.92)
bar(p,359,445,110,[(0,'0'),(.5,'0.04'),(1,'0.08')],'Speed / U')
text(p,420,431,'> 0.08');p.draw_line((460,433),(469,444),color=(0,0,0),width=.6)
p.draw_line((469,444),(463,441),color=(0,0,0),width=.6);p.draw_line((469,444),(468,438),color=(0,0,0),width=.6)
finish(d,'Fig_Pore_FlowField')

d=fitz.open();p=d.new_page(width=W,height=505)
text(p,8,16,'(a)');rect=(25,30,150,240.5263);field(p,'Device_speed',rect)
scale(p,'Device_speed',rect,10,'10 mm',.07,.98,below=True)
bar(p,20,274,164,[(0,'0'),(.5,'0.0225'),(1,'0.045')],'Axial speed (m/s)')
text(p,241,16,'(b)');chart(p,'Device_Profile_Chart',245,25)
text(p,8,319,'(c)')
chart(p,'Device_Loss_Chart',0,328)
finish(d,'Fig_Device_Flow')

src=json.loads((R/'data/fields/original_panel_provenance.json').read_text())['Fig_Tracer_Snapshots']
snap={f"{v['run']}_t{int(round(v['time_over_ideal']))}":v for v in src}
d=fitz.open();p=d.new_page(width=W,height=480)
for row,(stems,heading) in enumerate([(['E11_t1','E11_t3'],'G3, Pe = 3.6'),(['E12_t1','E12_t3'],'G3, Pe = 21.5'),(['E8_t2','E8_t10'],'G1, Pe = 138')]):
 y=17+row*132
 for col,stem in enumerate(stems):
  x=8+col*238;v=snap[stem]
  text(p,x,y+18,f"({chr(97+row*2+col)})")
  rect=(x,y+25,x+225,y+88.2);field(p,stem,rect);scale(p,stem,rect,20,'20 µm',.04,.98)
  residual=1-v['internal_average_from_report'];s=f'{residual:.4f}' if residual>=.0001 else '< 0.0001'
  text(p,x,y+106,'Internal residual: '+s)
bar(p,93,444,296,[(0,'0.001'),(1/3,'0.01'),(2/3,'0.1'),(1,'1')],'Unreplaced fraction (log scale)')
finish(d,'Fig_Tracer_Snapshots')
print('Composed three final-size Arial field plates.')
