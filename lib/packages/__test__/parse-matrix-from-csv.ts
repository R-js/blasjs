
/*function* generateLines(csv: string) {
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
            lineStart = i;
            (yield { text, line }) as void;
            line++
            continue;
        }
        i++;
    }
    if (lineStart < csv.length) {
        // one last return
        const text = csv.slice(lineStart, i);
        yield { text, line };
    }
}*/

function isSeparator(text: string, separator: string | RegExp): boolean {
    return (typeof separator === 'string') ? text?.[0] === separator : text[0].match(separator) !== null;
}

function absorbSeparator(text: string, i: number, separator: string | RegExp) {
    while (isSeparator(text[i], separator) && i < text.length) {
        i++;
    }
    return i; // absorbed till EOL
}

function isQuoted(text: string, i: number) {
    return text[i] === '"' || text[i] === '\'';
}

// text[i] is the position after the first quote
function absorbBetweenQuotes(text: string, i: number) {
    while (i < text.length
        &&
        (
            isQuoted(text, i) === false
            ||
            (isQuoted(text, i) && text[i - 1] === '\\')
        )
    ) {
        i++;
    }
    return i === text.length ? i : i + 1;
}

function absorbField(text: string, i: number, separator: string | RegExp, optionallyQuoted = true) {
    if (optionallyQuoted && isQuoted(text, i)) {// i am looking for a quote end
        return absorbBetweenQuotes(text, i + 1)
    }
    while (false === isSeparator(text[i], separator) && i < text.length) {
        i++;
    }
    return i;
}

function absorbComment(text: string, i: number) {
    if (text[i] === '#') {
        return text.length;
    }
    return i;
}

function absorbWhiteSpace(text: string, i: number) {
    while ('\t '.includes(text[i]) && i < text.length) {
        i++;
    }
    return i;
}

export type FieldType = 'field' | 'quoted-field' | 'comment' | 'column-name';

export function* generateFields(oneLine: { text: string, line: number }, separator: string | RegExp, optionallyQuoted = true) {
    let i = 0;
    const { text, line } = oneLine;
    let j = i;
    while (i < text.length) {
        i = absorbWhiteSpace(text, i); // absorb white space outside of fields
        j = absorbComment(text, i);
        if (j > i) {
            (yield { text: text.slice(i, j).trim(), line, s: i, e: j, type: 'comment' as FieldType }) as void;
            i = j;
        }
        j = absorbField(text, i, separator, optionallyQuoted);
        if (j > i) {
            (yield { text: text.slice(i, j).trim(), line, s: i, e: j, type: isQuoted(text, i) ? 'quoted-field' : 'field' as FieldType }) as void;
            i = j;
        }
        i = absorbSeparator(text, i, separator);
    }
}










