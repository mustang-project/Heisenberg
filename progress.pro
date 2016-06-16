pro progress, msg, cur, max
	writeu, -1, string(format='(%"\r",a,": ",i3,"%")', $
		msg, round(float(cur)/max*100))
	if cur eq max then print, '' ; newline
end
