default: clean

serve:
	jekyll serve

clean:

fresh:
	$(RM) -r _site .jekyll-cache
