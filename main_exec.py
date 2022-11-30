from datetime import datetime
import os
import subprocess
import sys
import time

def main():
	current_working_directory = str(os.getcwd())
	program = str(sys.argv[0])
	base_working_directory = os.path.dirname(os.path.abspath(program))
	log_filename = str(current_working_directory) + str('/log.txt')

	if os.path.exists(log_filename):
		append_write = 'a+' # append if already exists
	else:
		append_write = 'w' # make a new file if not

	file_object = open(log_filename, append_write)
	file_object.write(str("\n"))
	file_object.write(str("Log details for the run on: "))
	file_object.write(str(datetime.now()))
	file_object.write(str("\n"))
	file_object.close()
	
	code_name = "/main_code.py "

	global executable_string
	if (len(sys.argv)!=2):
		executable_string = str("python3 ") + str(base_working_directory) + code_name + str(current_working_directory) + str(" ") + str(base_working_directory) +  str(" | tee -a ") + str(log_filename)
	else:
		executable_string = str("python3 ") + str(base_working_directory) + code_name + str(current_working_directory) + str(" ") + str(base_working_directory) + str(" ") +  str(sys.argv[1]) + str(" | tee -a ") + str(log_filename)


if __name__ == "__main__":
	main()
	global executable_string
	os.system(executable_string)
	
