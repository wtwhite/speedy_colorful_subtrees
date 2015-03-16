%.csv: %.sql
	./select_to_csv.pl < $< > $@
