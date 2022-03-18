import pandas as pd
import re
def readscan(specfile,scanno):
    with open(specfile,'r') as fp:
            for i, line in enumerate(fp):
                if '#S {0:d}'.format(scanno) in line:
                    title = line
                    nlines = int(line.split()[-2])
                    nother = 0
                    start = False
                    for j, line  in enumerate(fp):
                        if '#L' in line:
                            cnames = re.split('\s{2,}',line[3:-1])   
                        elif '#' in line:
                            nother +=1
                            continue
                        elif (not '#' in line) \
                            and (not start):
                            headerindex = i+j-1
                            nstart = j
                            start = True
                        elif (not '#' in line) and  (
                            (start and '#C' in line) or 
                            (j == nlines+nother)):
                            nstop = j
                            break
                    break 
    nrows = nstop-nstart
    scan_info = {'columns':cnames[1:],'title':title}
    if nrows == 0:
        print('scan aborted')
        return([],scan_info)
    else:
        data = pd.read_csv(specfile,skiprows=headerindex+2,
                       nrows=nrows,delim_whitespace=True,header=None)                                    
        data.columns=cnames        
        return(data,scan_info)
               
def scanto(fd,flag):
	found = False
	while not found:
		t1 = fd.readline()
		if t1[0:len(flag)] == flag:
			found = True
	return t1

def scanthrough(fd,flag,stringin):
	found = True
	outstring = [stringin]
	while  found:
		t1 = fd.readline()
		if t1[0:len(flag)] != flag:
			found = False
		else:
			outstring.append(t1)
	return outstring

def getmnames(fd):
	mnames = scanthrough(fd,'#O',scanto(fd,'#O'))
	return mnames

def getmvals(fd,scanno,mnames):
	_ = scanto(fd,'#S {0:d}'.format(scanno))
	mvals = scanthrough(fd,'#P',scanto(fd,'#P'))
	marray = {}
	for ncol, vcol in zip(mnames,mvals):
		nvs = ncol.split()[1:]
		vvs = vcol.split()[1:]
		for tn, tv in zip(nvs,vvs):
			marray[tn] = tv
	return marray

def findXPCS(fd,flag):
	_  = scanto(fd,'#XPCS batch_name '+flag)
	lv = scanto(fd,'#S')
	vs = lv.split()
	return int(vs[1])-1