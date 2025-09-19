# Task # 1: To create a function that can convert the temperature ffrom celcius to fahrenheit

#temp_in_celcius = int(input("Enter the temperature in celcius: "))

# def convert_fahrenheit(temp):
#     converted_temp = temp * 9/5 + 32
#     return converted_temp

# converted = convert_fahrenheit(temp_in_celcius)

# print(f"The temperature from {temp_in_celcius}C is converted into {converted}F")


# Task # 2: To get maximum nuber from list of numbers

# list_of_numbers = [3, 5, 9, 92, 38, 92, 20]

# def max_number(numberlist):
#     largest = numberlist[0]
#     for x in numberlist:
#         if x > largest:
#             largest = x
#     return largest

# maximum_number = max_number(list_of_numbers)
    
# print(f"The greatest number in list is {maximum_number}")



# Task # 3

# user_input = input("Enter a sentence you want to add: ")

# def vowel_count(sentence):
#     count = 0
#     for character in sentence:
#         if character in 'aeiou':
#             count += 1
#     return count
    
# vowel = vowel_count(user_input)

# print(vowel)


# Task # 4 User credential is valid or not

existed_username = "iramjaved"
existed_password = "12345"


user_adding_username = input("Enter your username: ")
user_adding_password = input("Enter your password: ")


def validation_checking(username, password):
    istrue = False
    if existed_username == user_adding_username and existed_password == user_adding_password:
        istrue = True
    return istrue


valid = validation_checking(user_adding_username, user_adding_password)

if valid == True:
    print(f"Welcome: Username and password is correct.")
if valid != True:
    print(f"Sorry, your credentials are not valid.")









