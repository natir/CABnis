/*
Copyright (c) 2020 Pierre Marijon <pmarijon@mpi-inf.mpg.de>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
 */

/* crate use */
use thiserror::Error;

#[derive(Debug, Error)]
pub enum Error {
    #[error(
        "Reading of the file '{filename:}' impossible, does it exist and can be read by the user?"
    )]
    CantReadFile { filename: String },

    #[error("Creation/opening of the file '{filename:}' impossible, directory in path exist? can be write by the user?")]
    CantWriteFile { filename: String },

    #[allow(dead_code)]
    #[error("Error durring reading of file {filename:}")]
    ReadingError { filename: String },

    #[allow(dead_code)]
    #[error("Error durring writing of file {filename:}")]
    WritingError { filename: String },

    #[allow(dead_code)]
    #[error("If you get this error please contact the author with this message and command line you use: {name:?}")]
    NotReachableCode { name: String },
}
