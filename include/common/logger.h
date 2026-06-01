#ifndef LOGGER_H
#define LOGGER_H

/**
 *
 * @file logger.h
 *
 * @brief Prividema utility file for logging
 */

/**
 * @brief Represents the log level used in the program.
 *
 * This enumeration defines different flags that represent the log levels we can use while executing the program.
 * The user can choose the type of log by XORing the different flags in the enumeration.
 */
#include <stdint.h>
typedef enum {
	/// Detailed diagnostic information for developers.
	LOG_DEBUG = 0x1,

	/// General operational events showing that the system is functioning normally.
	LOG_INFO = 0x10,

	/// Something unexpected happened or might become a problem, but the system can continue running.
	LOG_WARN = 0x100,

	/// A failure occurred that prevented part of the system from functioning correctly.
	LOG_ERROR = 0x1000
} LogLevel;

/**
 * @brief Contains the log levels the user wants to display.
 *
 * This global variable contains the log levels the user wants to display.
 * The value is set by XORing the chosen flags in `LogLevel`.
 * Its default value is `0` which means no message will be displayed.
 */
extern uint64_t LOG_FLAG;

/**
 * @brief Sets the program's log flag.
 *
 * Updates the  global variable `LOG_FLAG`.
 * The value should be a XOR between chosen flags in `LogLevel`.
 *
 * @param flag The new value of the logFlag.
 */
void set_log_flag(unsigned char flag);

/**
 * @brief Logs a message with string formatting.
 *
 * @param level The log level.
 * @param fmt A va_list with the string to format and its arguments.
 *
 * @retval -1 if the level is `LOG_ERROR`.
 * @retval 0 otherwise.
 *
 * @note The formatting works like `printf`.
 */
int log_message(LogLevel level, const char* fmt, ...);

/**
 * @brief Logs a syscall error like perror but with `LOG_ERROR` flag.
 *
 * @param header The perror message header.
 *
 * @return Always returns `-1` as a syscall error is often fatal.
 */
int log_perror(const char* header);

/**
 * @brief Checks if a pointer is NULL. Log a syscall error like log_perror if so.
 *
 * @param ptr The pointer to test.
 *
 * @param header The perror message header.
 *
 * @retval -1 if and only if the pointer is NULL..
 * @retval 0 otherwise.
 */
int log_is_null(void* ptr, const char* header);

#endif  // LOGGER_H
