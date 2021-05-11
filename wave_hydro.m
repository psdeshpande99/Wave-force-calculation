disp('This is Group No.4')
disp('To calculate loads on cylinder:')
ip = input('Enter 1 for Ivagaki Method and 2 for Chakraborty method');
while (~(ip==1 || ip==2))
    ip = input('Please input valid numbers Enter 1 for Ivagaki Method and 2 for Chakraborty method');
end
if(ip==1)
        run('Wave_hydro_Ivagaki_uplusv_combine.m')
elseif(ip==2)
        run('Wave_hydro_Chakraborty_uandv_seperate.m')
end