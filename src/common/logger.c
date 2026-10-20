#include "logger.h"

#include <errno.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

uint8_t LOG_FLAG = LOG_ERROR | LOG_WARN | LOG_INFO | LOG_DEBUG;

// LOG_FLAG Setter
void set_log_flag(unsigned char flag) { LOG_FLAG = flag; }

// Get string representation of log level
const char* log_level_str(LogLevel level)
{
	switch (level & 0xff)
	{
		case LOG_DEBUG:
			return "DEBUG";
		case LOG_INFO:
			return "INFO ";
		case LOG_WARN:
			return "WARN ";
		case LOG_ERROR:
			return "ERROR";
		case LOG_USR1:
			return "USR1";
		case LOG_USR2:
			return "USR2";
		case LOG_USR3:
			return "USR3";
		case LOG_USR4:
			return "USR4";
	}
}

// The main logging function
int log_message(LogLevel level, const char* fmt, ...)
{
	// Skip the function if the level is not in LOG_FLAG
	if ((LOG_FLAG & level) == 0) return 0;

	// Get current time
	time_t t = time(NULL);
	struct tm tm_info;
	localtime_r(&t, &tm_info);  // thread-safe version

	char time_buf[20];
	strftime(time_buf, sizeof(time_buf), "%Y-%m-%d %H:%M:%S", &tm_info);

	// Print log level and timestamp
	printf("[%s] %s: ", time_buf, log_level_str(level));

	// Handle the variable arguments
	va_list args;
	va_start(args, fmt);
	vprintf(fmt, args);
	va_end(args);

	printf("\n");  // newline at the end

	return (level == LOG_ERROR) ? -1 : 0;
}

// Perror with logger
int log_perror(const char* header) { return log_message(LOG_ERROR, "%s : %s", header, strerror(errno)); }

int log_is_null(void* ptr, const char* header) { return (ptr == NULL) ? log_perror(header) : 0; }
