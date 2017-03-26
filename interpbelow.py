def interp(var,plist,target) :
    value = (var[1]-var[0])/(plist[1]-plist[0])*(target-plist[0])+var[0]
    return value