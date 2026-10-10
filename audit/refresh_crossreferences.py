from pathlib import Path
import re
ROOT=Path(__file__).resolve().parent.parent
for name,other,prefix in [('main','supplement','S-'),('supplement','main','M-')]:
    path=ROOT/(name+'.tex');s=path.read_text();aux=(ROOT/(other+'.aux')).read_text()
    labels={key:(num,page) for key,num,page in re.findall(r'\\newlabel\{([^}]+)\}\{\{([^}]*)\}\{([^}]*)\}',aux)}
    keys=sorted(set(re.findall(r'\\ref\*?\{'+re.escape(prefix)+r'([^}]+)\}',s)))
    for k in keys:
        if k not in labels:raise RuntimeError(f'Missing {other} label {k}')
    lines=['% BEGIN EMBEDDED CROSSREFERENCES','\\makeatletter']
    lines += [r'\@namedef{r@'+prefix+k+'}{{'+labels[k][0]+'}{'+labels[k][1]+'}{}{}{}}' for k in keys]
    lines += ['\\makeatother','% END EMBEDDED CROSSREFERENCES']
    s=re.sub(r'% BEGIN EMBEDDED CROSSREFERENCES.*?% END EMBEDDED CROSSREFERENCES','\n'.join(lines).replace('\\','\\\\'),s,flags=re.S)
    path.write_text(s)
    print(name, len(keys),'crossreferences refreshed')
