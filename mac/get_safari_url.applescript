#!/usr/bin/osascript

# About:
# A script that grabs all the URLs of every window open in Safari
# And saves them as html links in an .html file on the Desktop 
# (as well as a .txt file) 
# The file names are time-stamped

# Useage:
# ./get_safari_url.applescript

# Notes:
# borrowing from: 
# http://daringfireball.net/2003/02/save_and_restore_safari_urls
# http://www.macosxautomation.com/applescript/sbrt/sbrt-09.html
# http://www.cyberciti.biz/faq/mac-osx-applescript-run-shell-script/

# Get URLs of every open window in Safari:
tell application "Safari"
	set url_list to URL of every document
end tell

# convert url_list to text
set old_delim to AppleScript's text item delimiters

# set AppleScript's text item delimiters to newline:
set AppleScript's text item delimiters to "
"

set url_list to url_list as text
set AppleScript's text item delimiters to old_delim

# get timestamp and save to variable my_timestamp
set my_timestamp to time of (current date) as string

# save as text file (& concatenates strings)
set file_1 to (((path to desktop folder) as string) & "URL_list_" & my_timestamp & ".txt")
my write_to_file(url_list, file_1, false)

set file_2 to (((path to desktop folder) as string) & "URL_list_" & my_timestamp & ".html")

# get unix paths of files 
# (this is nec because applescript has its own stupid convention for paths with colons instead of slashes)
set file_1_path to POSIX path of file_1
set file_2_path to POSIX path of file_2

# convert text file to html links via shell awk 
do shell script "cat " & file_1_path & " | awk '{print \"<a href=\\\"\"$1\"\\\">\"$1\"</a><br>\"}' > " & file_2_path

# subroutine
on write_to_file(this_data, target_file, append_data)
	try
		set the target_file to the target_file as string
		set the open_target_file to open for access file target_file with write permission
		if append_data is false then set eof of the open_target_file to 0
		write this_data to the open_target_file starting at eof
		close access the open_target_file
		return true
	on error
		try
			close access file target_file
		end try
		return false
	end try
end write_to_file
