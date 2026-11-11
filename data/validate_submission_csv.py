# Faithful re-implementation of the site's client-side CSV validator (bundle fn Bv / tc / ec)
import sys
ALLOWED=set("ACDEFGHIKLMNPQRSTVWY")
def split_csv(line):
    out=[];cur='';q=False;i=0
    while i<len(line):
        c=line[i]
        if q:
            if c=='"':
                if i+1<len(line) and line[i+1]=='"': cur+='"'; i+=1
                else: q=False
            else: cur+=c
        else:
            if c=='"': q=True
            elif c==',': out.append(cur); cur=''
            else: cur+=c
        i+=1
    out.append(cur); return out
def validate(text, maxrows, mn, mx, req):
    import re
    text=text.lstrip('﻿')
    lines=re.split(r'\r?\n',text)
    nonblank=[l for l in lines if l.strip()]
    errs=[]
    if len(nonblank)>maxrows+1:
        return len(nonblank)-1,[(len(nonblank),'',f'最多提交 {maxrows} 条序列（含表头共 {maxrows+1} 行）')]
    hdr=[h.strip().lower() for h in split_csv(lines[0] if lines else '')]
    if 'name' not in hdr or 'sequence' not in hdr:
        return 0,[(1,'','表头必须包含 name 和 sequence')]
    oi,ui=hdr.index('name'),hdr.index('sequence')
    seen=set(); m=0
    for g in range(1,len(lines)):
        y=lines[g]
        if not y.strip(): continue
        j=split_csv(y)
        x=(j[oi] if oi<len(j) else '').strip()
        S=re.sub(r'\s+','',(j[ui] if ui<len(j) else '')).upper()
        p=g+1
        if not x and not S: continue
        m+=1
        if m>maxrows: return m,[(p,'',f'最多提交 {maxrows} 条序列（含表头共 {maxrows+1} 行）')]
        if not x or len(x)>64: errs.append((p,x,'名称须为 1 到 64 个字符')); continue
        if x in seen: errs.append((p,x,'名称重复')); continue
        seen.add(x)
        if not S: errs.append((p,x,'序列为空')); continue
        if any(ch not in ALLOWED for ch in S): errs.append((p,x,'含有非法氨基酸')); continue
        if len(S)<mn or len(S)>mx: errs.append((p,x,f'序列长度不符合要求：{req}'))
    if m==0: errs.append((1,'','没有序列'))
    return m,errs
if __name__=='__main__':
    txt=open(sys.argv[1],encoding='utf-8-sig').read()
    n,e=validate(txt,20,40,100,'minibinder 为 40—100 aa')
    print('file:',sys.argv[1],'| parsed sequences:',n,'| errors:',e if e else 'NONE -> would upload')
