
function* generateLines(csv: string) {
    // parse for line endings
    let lineStart = 0;
    let line = 1;
    let i = 0;
    while (i < csv.length) {
        // if we see a \r\n treat it as one
        if (csv[i] === '\r') {
            const text = csv.slice(lineStart, i);
            (yield { text, line }) as void;
            line++;
            if (csv[i + 1] === '\n') {
                i += 2;
            }
            else {
                i += 1;
            }
            lineStart = i;
            continue;
        }
        if (csv[i] === '\n') {
            const text = csv.slice(lineStart, i);
            i++;
            line++;
            lineStart = i;
            (yield { text, line }) as void;
        }
    }
    if (lineStart < csv.length){
        // one last return
        const text = csv.slice(lineStart, i);
        yield { text, line};
    }
}

function* generateFields(oneLine: {text: string, line: number}, seperator: string| RegExp ){
    let col = 1;


}



