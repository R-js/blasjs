
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
    if (lineStart < csv.length) {
        // one last return
        const text = csv.slice(lineStart, i);
        yield { text, line };
    }
}

function isSeparator(text: string, separator: string | RegExp): boolean {
    return (typeof separator === 'string') ? text[0] === separator : text[0].match(separator) !== null;
}

function absorbSeparator(text: string, i: number, separator: string | RegExp) {
    while (isSeparator(text[i], separator) && i < text.length) {
        i++;
    }
    return i; // absorbed till EOL
}

function absorbField(text: string, i: number, separator: string | RegExp) {
    while (false === isSeparator(text[i], separator) && i < text.length) {
        i++;
    }
    return i;
}

function* generateFields(oneLine: { text: string, line: number }, separator: string | RegExp) {
    let i = 0;
    const { text, line } = oneLine;
    while (i < text.length) {
        i = absorbSeparator(text, i, separator);
        let j = absorbField(text, i, separator);
        if (j > i){
            (yield { field: text.slice(i,j), line, i}) as void;
            i = j;
        }
    }
}

export default function parse(text: string, hasColumnNames = true, hasRowNames = true){


}



