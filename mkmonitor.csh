#!/bin/sh
cat << END
#!/bin/sh
if [ \$# -gt 0 ]
then
	option=\$1
else 
	option='on'

fi
if [ \$option != "on" ] 
then
	touch monitor.stop
	exit
fi
rm -f *.monitor
$PWD/monitor.tk &
exit
END
exit
