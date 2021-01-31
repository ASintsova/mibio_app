import streamlit as st
import pandas as pd


def load_files(uploaded_files):
    df_list = []
    for file in uploaded_files:
        file.seek(0)
        df_list.append(pd.read_csv(file, index_col = 0))
    return pd.concat(df_list)

def main():
    st.title("TNSeq Data Exploration and Analysis")

    # Load the data

    files = list(st.file_uploader(label="Add count files", accept_multiple_files=True))
    if not files:
        st.stop()
    df = load_files(files)
    st.write(df.sample(10))

if __name__ == "__main__":
    main()