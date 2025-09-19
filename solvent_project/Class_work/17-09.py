# Task # 1: 

# def username_sanitized(username):
#     usernamee = username.lower()
#     usernamee = username.replace(" ", "")
#     sanitized = ""
#     for char in username:
#         if char.isalnum():
#             sanitized += char
#     return sanitized

# def get_username():
#     username_by_user = input("Enter your username: ")
#     return username_by_user

# input_by_user = get_username()

# sanitized_username = username_sanitized(input_by_user)


# print(f"The sanitized username is {sanitized_username}")


# Task # 2 

# def count_word(text):
#     wordscount = text.split()
#     return len(wordscount)

# def get_text():
#     text_by_user = input("Enter text: ")
#     return text_by_user

# user_input = get_text()

# word_count_of_user_text = count_word(user_input)

# print(f"The word count of text you provided is {word_count_of_user_text}")


# Task # 3

# def email_checking(email):
#     istrue = False
#     if email.endswith("@email.com") or email.endswith("@gmail.com") or email.endswith("@yahoo.com"):
#         istrue = True
#     return istrue

# def get_email():
#     email_by_user = input("Enter email: ")
#     return email_by_user

# user_email = get_email()

# check_email = email_checking(user_email)

# if check_email == True:
#     print(f"Passed")
# if check_email != True:
#     print(f"Your email is not in right format. Write in format @email.com")


# Task # 4 

def remove_words_from_text(given_text):
    new_text = given_text.remove()
    return new_text

def get_data_text():
    data_by_user = input("Enter line: ")
    return data_by_user

get_data = get_data_text()

removed_words = remove_words_from_text(get_data)

print(f"Finalised line: {removed_words}")



