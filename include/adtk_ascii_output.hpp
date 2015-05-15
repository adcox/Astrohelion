#ifndef __H_ASCII_
#define __H_ASCII_

/* 
 *	ASCII Escape Sequences for text coloring. Usage:
 *		printf( RED "Here is some red text!" RESET);
 */

#define RESET   "\033[0m"					/**< Reset ASCII text to default values */
#define BLACK   "\033[30m"					/**< Make ASCII text black */
#define RED     "\033[31m"					/**< Make ASCII text red */
#define GREEN   "\033[32m"					/**< Make ASCII text green */
#define YELLOW  "\033[33m"					/**< Make ASCII text yellow */
#define BLUE    "\033[34m"					/**< Make ASCII text blue */
#define MAGENTA "\033[35m"					/**< Make ASCII text magenta */
#define CYAN    "\033[36m"					/**< Make ASCII text cyan */
#define WHITE   "\033[37m"					/**< Make ASCII text white */
#define BOLDBLACK   "\033[1m\033[30m"		/**< Make ASCII text bold and black */
#define BOLDRED     "\033[1m\033[31m"		/**< Make ASCII text bold and red */
#define BOLDGREEN   "\033[1m\033[32m"		/**< Make ASCII text bold and green */
#define BOLDYELLOW  "\033[1m\033[33m"		/**< Make ASCII text bold and yellow */
#define BOLDBLUE    "\033[1m\033[34m"		/**< Make ASCII text bold and blue */
#define BOLDMAGENTA "\033[1m\033[35m"		/**< Make ASCII text bold and magenta */
#define BOLDCYAN    "\033[1m\033[36m"		/**< Make ASCII text bold and cyan */
#define BOLDWHITE   "\033[1m\033[37m"		/**< Make ASCII text bold and white */

#endif
